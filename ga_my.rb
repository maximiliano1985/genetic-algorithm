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
    def initialize( i_o = {} )
        raise ArgumentError, "Need an Integer instead of #{dim.class}" unless dim.kind_of? Integer
        raise ArgumentError, "Need an Hash instead of #{i_o.class}" unless i_o.kind_of? Hash
        @dimension = i_o.size # is number of chromosomes in the population
        max_v = 0.0
        # find the greatest number in the interval of values  the first population
        i_o.each_value { |v| max_v = v.max if v.max > max_v}
        num_i = max_v.to_i # the integer part
        num_d = ( (max_v-max_v.to_i)*1000 ).to_i# the decimal part
        @nbit = [ num_i.to_s(2).size , num_d.to_s(2).size ].max # set the number of bit necessary to encode in a binary gene
        @size = i_o.size     # is the problem size, the domain dimention of objective function
        @population = [] # initialize the population, is an hash which keys are: :chromosome and :fitness
        rr = Random.new()
        dim.times do # for each chromosome of the initial population do...
            cr = []
            # generate the gene randomly (it must lie in the domain defined by i_o)
            i_o.each_value{ |v| cr << rr.rand(v.min..v.max) } # generates a random number from a uniform distribution
            @population << { :chromosome => cr } 
        end
    end
    
    # Is the difference of the function value for the best and the worst chromosome 
    # Input: array of hashes. Output: float
    def norm
        pop = @population.sort{ |x, y| y[:fitness] <=> x[:fitness] }
        best = pop.first
        worst = pop.last
        return ( best[:fitness] - worst[:fitness] ).abs
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
          :tol        => 0.001, # accurancy of the solution 
          :p_mutation => 0.2,   # probability of mutation
          :p_crossover=> 0.8,   # probability of cross over
          :i_o        => {},    # inverval of values used for define the first population
          :npop       => 50,    # number of population to be computed
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
     
    def loop( ary = nil )
        raise ArgumentError, "Block needed" unless block_given?
        # evaluates the cromosomes and converts them into a string of bits
        @population.each do |c|
            c[:fitness]   = yield( c[:chromosome] ) # evaluates the chromosome
            c[:bitstring] = encode( c[:chromosome] ) # converts the chromosome into a string of bits
        end
        
        until converged? or @iteration > @cfg[:npop]
            selected = Array.new( @population.size ){ |i| selection( @population ) }
            bit_selected = []
            selected.each{ |v| bit_selected << v[:bitstring] }
            # childs_a is an array of floats. @population.size is used to set the number of chromosomes in the new generation
            childs = evolve( bit_selected, @population.size, @cfg[:p_crossover], @cfg[:p_mutation] ) 
            # child is converted into an array of hashes, each with keys: :chromosome, :bitstring , :fitness
            @population = [] # reset the population and then update it
            childs.each do |c|
                @population << {
                 :bitstring  => c,    
                 :chromosome => decode( c ), # converts the string of bits into an array of floats
                 :fitness    => yield( c.decode )# evaluates the array of floats  
                }
            end # childs do

            puts "#{@iteration}st generation, best: #{best[:chromosome]} ---> #{best[:fitness]}"
            
            # these lines ar used to do a convergence plot, i.e. all the fitnesses for the current population
            if @cfg[:pconv] == true
                sdata = [] if @iteration == 0 # at first iteration initialize the matrix containing simplex data
                @population.each{ |v| sdata << v[:fitness] }
                sleep 0.1
                # a. initialize the data sets for the plot
                @gp.new_series
                # b. fill the data sets
                sdata.each_with_index { |v,i| @gp.series[ i.to_s.to_sym ] << [ @iteration , v ] }
                # c. close the data sets
                @gp.series.close
                # d. plot the data sets
                @iteration == 0 ? @gp.plot , "with points"    :    @gp.replot , "with lines"
            end # if @cfg
            @iteration += 1
        end # until converged
        return @population.sort!{ |x, y| y[:fitness] <=> x[:fitness] }.first
    end # def loop
    
    private
    # Input: @cfg. Output: boolean value
    def converged?
      n = @population.norm
      if n then
        @population.norm < @cfg[:tol] # this returns true or false
      else
        false
      end
    end
    
    # Input: array of bit strings. Output: array of bit strings
    def evolve(selected, pop_size, p_cross, p_mut)
        children = []
        selected.each_with_index do |p1, i|
            # crossing over
            p2 = (i.modulo(2) == 0) ? selected[i+1] : selected[i-1] # i.modulo(2) is the reminder of the division i / 2 
            p2 = selected[0] if i == selected.size - 1
            ch1 = {} ; ch2 = {} ; ch3 = {}
            ch1[:chromosome] , ch2[:chromosome] = crossover( p1[:chromosome], p2[:chromosome], p_cross )
            children << ch1
            childern << ch2
            
            # mutation
            p3 = (i.modulo(2) == 1) ? selected[i+1] : selected[i-1] # p3 is a chromosome not yet used
            ch3[:chromosome] = mutation( p3[:chromosome], p_mut )
            children << ch3
            break if children.size >= pop_size
        end
        return children
    end
    
    # compares two chromosomes and selects the one with the best fitting.
    # This is a binary tournament.
    # Input: array of hashes. Output: array of hashes
    def selection(pop)
        selected = []
        raise "Error in selection: input must be a Array instead of a #{pop.class}" unless pop.kind_of? Array
        i , j = rand( pop.size ) , rand(pop.size)
        j = rand(pop.size) while j == i # if unfortunatly j = i, evaluates j again
        selected << (pop[i][:fitness] > pop[j][:fitness]) ? pop[i] : pop[j]
        return selected
        # the size of selected is the half of the population size
    end
    
    # the chromosome is a string of '0' and '1', rate [0,1]. mutant is also a string
    # Input: a bit string, Output: a bit string
    def mutation( bitstring , rate )
        raise "Error in mutation: input must be a String instead of a #{bitstring.class}" unless bitstring.kind_of? String
        mutant = ""
        bitstring.size.times do |i|
            gene = bitstring[i]
            # change the bit value only if the rand [0,1] is minor than the mutation probability
            mutant << ((rand < rate) ? ((gene == '1') ? "0" : "1") : gene)
        end # chromosome
        return mutant
    end # def mutation
    
    # both father and mather are strings, rate [0,1] is the crossover probability
    # the crossover returns two childs.
    # Input: 2 bit strings + 1 float. Output: 2 bit strings
    def crossover(father, mother, rate)
        raise "Error in crossover: input must be a String instead of a #{father.class}" unless father.kind_of? String
        return father , mother if rand >= rate # don't do the cross over if rand is maior or equal than the crossover probability 
        point = = (rand*mother.size).to_i # sets the crossover point randomly
        return father[0..point] + mother[point..(mother.size)] , mother[0..point] + father[point..(father.size)]
    end
    
    # this is used to convert the input values into a binary string (the chromosome)
    # Input: array of floats. Output: one bit string
    def encode( ary = [] )
        ary_b = []
        ary.each { |v|
            ary_b << v.to_i.to_bin( @nbit )
            ary_b << ( ( (v-v.to_i)*1000 ).to_i ).to_bin( @nbit )
        }
        return ary_b*""
    end
    
    # this is used to decode the binary chromosome string  into an array of floats
    # Input: one bit string. Output: array of floats
    def decode( str )
        ng = str.size / @nbit # is the number of genes in one chromosome
        raise "Error in the chromosome decodification: you have #{ng} genes of #{@nbit} bits, and #{str.size} bits in the chromosome" unless ng * @nbit == str.size
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
  end # class Optimizer
end # module GA

if __FILE__ == $0 then
  # Test function
  f = lambda {|p| p[0]**2 + 3*p[1]**2+10} # a trivial parabola
  #f = lambda { |p| p[0] ** 2 - 4 * p[0] + p[1] ** 2  -p[1] - p[0] * p[1]}
  #f = lambda { |p| ( 1 - p[0] ) ** 2 + 100 * ( p[1] - p[0] ) ** 2 } # Rosenbroke function
  
  # Instantiate the optimizer, with tolerance and dimension
  pop = GA::Population.new( dom={:X =>[5,10], :Y=>[-10.23,54.234]}  )  
  p pop.population
  pop_b = []
  pop.population.each{ |v| pop_b << pop.encode(v[:chromosome]) } 
  p pop_b
  #opt = NMM::Optimizer.new( :dim => 3, :tol => 1E-5,:niter => 50, :exp_f => 2,:cnt_f => 0.5)
  
  # Define the starting points, i.e. the first guess simplex
  #opt.start_points = [
  #  Vector[0.1,0.1],
  #  Vector[1.2,0.0],
  #  Vector[0.0,0.8] 
  #]
  #opt.loop {|p| f.call(p)}

end
