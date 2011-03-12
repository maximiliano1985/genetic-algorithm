#!/usr/bin/env ruby
#  ga_my
#
#  Created by Carlos Maximiliano Giorgio Bort on 2011-03-06.
#  Copyright (c) 2011 University of Trento. All rights reserved.
#

require 'rubygems'
require 'gnuplotr'

class Vector
  if RUBY_VERSION.split(".").join.to_i < 190
    # Add method for division under Ruby 1.8.x
    def /(other); self * (1.0 / other); end
  end
  # More compact inspect version
  def inspect; "V[#{self.to_a * ","}]"; end
end
class Integer
    def to_bin(size) # size is the number of bits used 
        num = self.to_i # the integer part
        bin = num.to_s(2)
        bin.size <= size ? inc = (size - bin.size) : (raise "Use more bits for a proper binary convertion")
        inc.times { bin = "0" +  bin }
        return bin
    end
end

module GA
  class Population
    attr_reader :population

    # i_o is the inverval of values used for define the first population
    # values for the parameters of the first population
    def initialize(i_o = {})
        raise ArgumentError, "Need an Integer instead of #{dim.class}" unless dim.kind_of? Integer
        raise ArgumentError, "Need an Hash instead of #{i_o.class}" unless i_o.kind_of? Hash
        @dimension = i_o.size # is number of cromosomes in the population
        max_v = 0.0
        # find the greatest number in the interval of values  the first population
        i_o.each_value { |v| max_v = v.max if v.max > max_v}
        num_i = max_v.to_i # the integer part
        num_d = ( (max_v-max_v.to_i)*1000 ).to_i# the decimal part
        @nbit = [ num_i.to_s(2).size , num_d.to_s(2).size ].max # this is used to set the number of bit necessary to encode in a binary gene the inputs
        @size = i_o.size     # is the problem size, the domain dimention of objective function
        @population = [] # initialize the population, is an hash which keys are: :cromosome and :fitness
        rr = Random.new()
        dim.times do # for each cromosome of the initial population do...
            cr = []
            # generate the gene randomly (it must lie in the domain defined by i_o)
            i_o.each_value{ |v| cr << rr.rand(v.min..v.max) } # generates a random number from a uniform distribution
            @population << { :cromosome => cr } 
        end
        @ready = false # tell then the simplex is ready for the optimization
    end
    
    # Performs the anslyis of the population. It sorts the cromosomes in ascending order
    # acordingly to theirs fitness.
    def reorder
        return if @ready
        @population = @population.sort {|a,b| a[:fitness]<=>b[:fitness] }
    end
    
    # this is used to convert the input values into a binary string (the cromosome)
    # input: Array, output: String
    def encode(ary=[])
        ary_b = []
        ary.each { |v|
            ary_b << v.to_i.to_bin( @nbit )
            ary_b << ( ( (v-v.to_i)*1000 ).to_i ).to_bin( @nbit )
        }
        return ary_b*""
    end
    
    # this is used to decode the binary cromosome string  into an array of floats
    # input: String, output: Array
    def decode(str)
        ng = str.size / @nbit # is the number of genes in one cromosome
        raise "Error in the cromosome decodification: you have #{ng} genes of #{@nbit} bits, and #{str.size} bits in the cromosome" unless ng * @nbit == str.size
        cr = []
        dots = "" 
        @nbit.times{ dots += "." } # generates a string with @nbit dots: "....."
        dots = "(" + dots + ")"    # adds the parentesys: "(.....)"
        d = Regexp.new( dots )     # converst dots into a regulare expression: /(.....)/
        str_a = str.split( d )     # splits the string: eg. "00000111112222233333" -> ["", "00000", "", "11111", "", "22222", "", "33333"]
        str_a = str_a.delete("")   # ["00000", "11111", "22222", "33333"]
        str_a.each_index{ |i| cr << ( str_a[2*i]+"."+str_a[2*i+1] ).to_i(2) if i <= str_a.size-2} # this returns: [ 00000.11111 , 22222.33333 ], each element in the array is in decimal codification
        return cr
    end
    
    def performance
    end
         
    end
  end # class Population

  # Class that implements a general n-dimensional Nelder-Meade Method (NMM).
  # @author Paolo Bosetti
  class Optimizer
    attr_reader :simplex, :status, :iteration
    # The +Optimizer+ gets initialized with the following defaults:
    #     :dim   => 2,
    #     :exp_f => 1.5, 
    #     :cnt_f => 0.5,
    #     :tol   => 0.001
    # If no +args+ is specified, these are the defaults. Otherwise, the keys that 
    # are passed are merged with those defaults.
    # @param [Hash] args a +Hash+ of initialization values
    def initialize(args = {})
        @cfg = {
          :tol        => 0.001, # the accurancy of the solution 
          :p_mutation => 0.2,   # the probability of mutation
          :p_crossover=> 0.8,   # the probability of cross over
          :i_o        => {},
          :npop       => 50,     # the number of population to be computed
          :pconv      => true,
          :plotopt    => {:title  => 'Genetic Algorithm Convergence',
                          :xlabel => 'No. iteration',
                          :ylabel => 'Objective function value',
                          :yrange => [ -10 , 10 ],
                          :grid   => "set grid"
                         }
        }
        raise "Error with the assigned mutation probability:\n it is #{@cfg[:p_mutation]} but must be 0 <= p_mutation <= 1 " unless @cfg[:p_mutation] >= 0 and @cfg[:p_mutation] <= 1 
        raise "Error with the assigned crossover probability:\n it is #{@cfg[:p_crossover]} but must be 0 <= p_crossover <= 1 " unless @cfg[:p_crossover] >= 0 and @cfg[:p_crossover] <= 1 
        @cfg.merge! args
        @population = Population.new(@cfg[:i_o])
        @start_points = []
        @status = :filling
        @iteration = 0
        if @cfg[:pconv] == true # this is the plot
            @gp = GNUPlotr.new
            # enable command history recording
            @gp.record = true
            # Issue raw gnuplot commands
            @gp.raw @cfg[:plotopt][:grid]
            # Some magic mapping works too:
            @gp.set_grid
            @gp.set_title @cfg[:plotopt][:title]  , :font => "Times New Roman,18"
            @gp.set_xlabel @cfg[:plotopt][:xlabel], :font => "Times New Roman,18"
            @gp.set_ylabel @cfg[:plotopt][:ylabel], :font => "Times New Roman,18"
            @gp.set_xrange( 0 .. @cfg[:niter])
            @gp.set_yrange(@cfg[:plotopt][:yrange][0] .. @cfg[:plotopt][:yrange][1])
        end # if @cfg
    end
     
    def evolve
    end
    
    # compares two cromosomes and selects the one with the best fitting.
    # This is binary tournament
    def selection(pop)
        i , j = rand( pop.size ) , rand(pop.size)
        j = rand(pop.size) while j == i # if unfortunatly j = i, evaluates j again
        return (pop[i][:fitness] > pop[j][:fitness]) ? pop[i] : pop[j]
    end
    
    # the cromosome is a string of '0' and '1', rate [0,1]. mutant is also a string
    def mutation(cromosome, rate)
        mutant = ""
        cromosome.size.times do |i|
            gene = cromosome[i]
            # change the bit value only if the rand [0,1] is minor than the mutation probability
            mutant << ((rand < rate) ? ((gene == '1') ? "0" : "1") : gene)
        end # cromosome
        return mutant
    end # def mutation
    
    # both father and mather are strings, rate [0,1] is the crossover probability
    # the crossover returns two childs 
    def crossover(father, mother, rate)
        
    end
    
    def invertion
    end
  end # class Optimizer
end # module GA

if __FILE__ == $0 then
  pop = GA::Population.new( 10 , dom={:A =>[0,1], :B=>[-10.23,54.234]}  )  
  p pop.population
  pop_b = []
  pop.population.each{ |v| pop_b << pop.encode(v[:cromosome]) } 
  p pop_b
    
  # Test function
  #f = lambda {|p| p[0]**2 + 3*p[1]**2+10} # a trivial parabola
  #f = lambda { |p| p[0] ** 2 - 4 * p[0] + p[1] ** 2  -p[1] - p[0] * p[1]} # a non trivial parabola
  #f = lambda { |p| ( 1 - p[0] ) ** 2 + 100 * ( p[1] - p[0] ) ** 2 } # the non trivial Rosenbroke function
  # Instantiate the optimizer, with tolerance and dimension (it is the dimension
  # of the simplex, so number of parameters + 1 !!!)
  #opt = NMM::Optimizer.new( :dim => 3, :tol => 1E-5,:niter => 50, :exp_f => 2,:cnt_f => 0.5)
  
  # Define the starting points, i.e. the first guess simplex
  #opt.start_points = [
  #  Vector[0.1,0.1],
  #  Vector[1.2,0.0],
  #  Vector[0.0,0.8] 
  #]
  #opt.loop {|p| f.call(p)}

end
