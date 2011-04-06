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
    def to_bin(siz) # size is the number of bits used 
        bin = self.to_s(2)
        bin.size <= siz ? inc = (siz - bin.size) : (raise "Use more bits for a proper binary convertion, you used #{size} bits for a #{bin.size} binary number")
        inc.times { bin = "0" +  bin }
        return bin
    end
end

module GA
  class Population
    attr_reader :population , :nbit

    # i_o is the inverval of values used for define the first population
    # values for the parameters of the first population
    def initialize( dim, i_o = {}, prec = 1E-3 )
        raise ArgumentError, "Need an Hash instead of #{i_o.class}" unless i_o.kind_of? Hash
        @population = [] # initialize the population, is an hash which keys are: :chromosome and :fitness
        dim.times do # for each chromosome of the initial population do...
            cr = []
            # generate the gene randomly (it must lie in the domain defined by i_o)
            i_o.each_value{ |v| cr << rand()*(v.max-v.min) + v.min }
            @population << { :chromosome => cr }
        end
        @prec = prec # the precision, i.e. the number considered at the left of the comma
    end
    
  end # class Population

  # Class that implements a general n-dimensional Nelder-Meade Method (NMM).
  # @author Paolo Bosetti
  class Optimizer
    attr_reader :simplex, :status, :iteration
    def initialize(args = {})
        @cfg = {
          :toll       => 1E-3 , # accurancy of the solution 
          :nbit       => 10   , # number of bits used to encode the cromosomes into a binary string
          :p_mutation => 0.2  , # probability of mutation
          :p_crossover=> 0.8  , # probability of cross over
          :i_o        => {}   , # inverval of values used for define the first population
          :npop       => 50   , # number of population to be computed
          :ncr        => 100  , # number of chromosomes in each population
          :pconv      => true ,
          :nelitist   => 3    , # the 'n' best chromosomes that will automatically be copied in the new population
          :plotopt    => {:xlabel => 'No. iteration',
                          :ylabel => 'Objective function value',
                          :yrange => [ -10 , 10 ],
                          :grid   => "set grid"
                         }
        }
        @cfg[:plotopt][:title] = "Population size: #{@cfg[:ncr]} chromosomes"
        raise "Error with the assigned mutation probability:\n it is #{@cfg[:p_mutation]} but must be 0 <= p_mutation <= 1 " unless @cfg[:p_mutation] >= 0 and @cfg[:p_mutation] <= 1 
        raise "Error with the assigned crossover probability:\n it is #{@cfg[:p_crossover]} but must be 0 <= p_crossover <= 1 " unless @cfg[:p_crossover] >= 0 and @cfg[:p_crossover] <= 1 
        @cfg.merge! args
        @nbit = @cfg[:nbit] #@pop.max_bit # is the number f bits required to encode into a binary string the chromosome
        @pop = Population.new( @cfg[:ncr] , @cfg[:i_o])
        @population = @pop.population
        @start_points = []
        @status = :filling
        @iteration = 0
        @best = []
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
            @gp.set_xrange( 0 .. @cfg[:npop])
            #@gp.set_yrange(@cfg[:plotopt][:yrange][0] .. @cfg[:plotopt][:yrange][1])
        end # if @cfg
    end 
     
    def loop( ary = nil )
        raise ArgumentError, "Block needed" unless block_given?
        # evaluates the cromosomes and converts them into a string of bits
        @population.each do |c|
            c[:fitness]   = yield( c[:chromosome] )  unless  c[:fitness]   # evaluates the chromosome
            c[:bitstring] = encode( c[:chromosome] ) unless c[:bitstring] # converts the chromosome into a string of bits
        end
        
        until converged? or @iteration > @cfg[:npop]
            @sorted = @population.sort!{ |x, y| x[:fitness] <=> y[:fitness] }
            selected = selection( @population )
            bit_selected = []
            selected.each{ |v| bit_selected << v[:bitstring] }
            # @population.size is used to set the number of chromosomes in the new generation
            childs = evolve( bit_selected, @population.size-@cfg[:nelitist] , @cfg[:p_crossover], @cfg[:p_mutation] ) 
            # child is converted into an array of hashes, each with keys: :chromosome, :bitstring , :fitness
            @population = @sorted[ 0 .. @cfg[:nelitist]-1 ] # reset the population and then update it
            childs.each do |c|
                dec =  decode( c[:bitstring] )
                @population << {
                    :bitstring  => c[:bitstring],
                    :chromosome => dec,           # converts the string of bits into an array of floats   
                    :fitness    => yield( dec )   # evaluates the array of floats  
                }
            end # childs do
            @sorted = @population.sort!{ |x, y| x[:fitness] <=> y[:fitness] }
            @best << @sorted.first
            puts "#{@iteration}th generation, best: #{@best.last.inspect}" ##########
            puts "#{@iteration}th generation, worst: #{ @sorted.last.inspect }" ###########
            "Maximum number of iteration reached: #{@cfg[:npop]}" if @iteration == @cfg[:npop] 
            puts "_________________________________________________________"
            
            # these lines ar used to do a convergence plot, i.e. all the fitnesses for the current population
            if @cfg[:pconv] == true
                
                # a. initialize the data sets for the plot
                @gp.new_series(:population)
                
                # b. fill the data sets
                if @iteration == 0 # initialize the matrix containing simplex data
                    sdata = []
                    it = []
                end
                it << @iteration # array of integers

                ### these lines are usefull to plot the evolution of entire population
                #f_a = []
                #@population.each{ |v| f_a << v[:fitness] } # array of array of floats
                #sdata << f_a
                # b. fill the data sets
                #it.each do |i|
                #    #sdata[i].each do |v| 
                #    sdata.each do |v| 
                #        @gp.series[:population] << [ i , v ]
                #    end
                #end
                
                ### these lines are used to track the best chromosome in the population
                sdata << @best[-1][:fitness]
                it.each do |i|
                    @gp.series[:population] << [ i , sdata[i] ]
                end
                # c. close the data sets
                @gp.series[:population].close
                # d. plot the data sets
                if @iteration == 0
                    @gp.plot :population   ,  "with points lt 9 pt 2 notitle"
                    @gp.plot :population   ,  "with line lt 9 pt 2 notitle"
                else
                    @gp.replot :population , "with points lt 9 pt 2 notitle"
                    @gp.plot :population   ,  "with line lt 9 pt 2 notitle"
                end
            end # if @cfg
            @iteration += 1
        end # until converged
        return @best[-1]
    end # def loop
    
    private
    
    # The solution converges if the fitness for the best chromosome of the latter 3 population is the same
    # Input: array of hashes. Output: boolean value
    def converged?
      if @iteration > 1 
          xx = 0.0
          3.times{ |p| xx += @sorted[p-1][:fitness]**2 }
          if xx**0.5 <= @cfg[:toll]
              true
          else
              false
          end # if xx
      end # if @iteration
    end # converged
    
    # Input: array of bit strings. Output: array of bit strings
    def evolve( selected, pop_size, p_cross, p_mut ) 
        children = []
        selected.each_with_index do |p1, i|
            # crossing over
            ############################## non me piase sto metodo di scelgiere moglie e marito
            # 1 ::: p2 = (i.modulo(2) == 0) ? selected[i+1] : selected[i-1] # i.modulo(2) is the reminder of the division i / 2 
            #       p2 = selected[0] if i == selected.size - 1
            # 2 ::: i == selected.size-1 ? p2 = selected[0] : p2 = selected[i+1] 
            # 3 ::: random choise:
            j = rand(selected.size)
            j = rand(selected.size) while j == i
            p2 = selected[j]
            ch1 = {} ; ch2 = {} ; ch3 = {}
            ch1[:bitstring] , ch2[:bitstring] = crossover( p1, p2, p_cross )
            children.concat( [ ch1 , ch2 ] )
            
            # mutation
            #i.modulo(2) == 1 ? p3 = selected[i+1] : p3 = selected[i-1] # p3 is a chromosome not yet used
            k = rand(selected.size)
            k = rand(selected.size) while k == j and k == i # random choise of the chromosome that might mutate
            p3 = selected[k]
            ch3[:bitstring] = mutation( p3, p_mut )
            children << ch3
            break if children.size >= pop_size
        end
        return children
    end
    
    # compares two chromosomes and selects the one with the best fitting.
    # the size of selected is the half of the population size
    # This is a binary tournament.
    # Input: array of hashes. Output: array of hashes
    def selection(pop)
        raise "Error in selection: input must be a Array instead of a #{pop.class}" unless pop.kind_of? Array
        sel = []
        pop.size.times do
            i , j = rand( pop.size ) , rand( pop.size )
            j = rand( pop.size ) while j == i # if unfortunatly j = i, evaluates j again
            (pop[i][:fitness] < pop[j][:fitness]) ? sel << pop[i] : sel << pop[j]
        end
        return sel
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
    def crossover( father, mother, rate )
        raise "Error in crossover: father must be a String instead of a #{father.class}" unless father.kind_of? String
        raise "Error in crossover: mother must be a String instead of a #{mother.class}" unless mother.kind_of? String
        # don't do the cross over if rand is maior or equal than the crossover probability 
        return father , mother if rand >= rate
        raise "Error in crossover, father and mother must have the same dimension" unless father.size == mother.size
        point = 0.0
        point = rand(mother.size) while point == 0 or point == mother.size # sets the crossover point randomly
        return father[0..point-1] + mother[point..(mother.size)] , mother[0..point-1] + father[point..(father.size)] 
    end
    
    # this is used to convert the input values into a binary string (the chromosome)
    # Input: array of floats. Output: one bit string
    def encode( ary = [] )
        ary_b = ""
        ary.each { |v|
            i_b =  v.to_i.abs.to_bin( @nbit )
            # how to encode the sign: the first bit is 0 if i_b is >= 0, and 1 if it's < 0
            v.to_i >= 0 ? i_b[0] = "0" : i_b[0] = "1"
            d_b = ( ( (v-v.to_i)/@cfg[:tol] ).to_i.abs ).to_bin( @nbit )
            ary_b += i_b+d_b
        }
        return ary_b
    end
    
    # this is used to decode the binary chromosome string  into an array of floats
    # Input: one bit string. Output: array of floats
    def decode( str )
        ng = str.size / @nbit # is the number of genes in one chromosome
        raise "Error in the chromosome decodification: you have #{ng} genes of #{@nbit} bits, and #{str.size} bits in the chromosome" unless ng * @nbit  == str.size
        cr = []
        dots = "" 
        @nbit.times{ dots += "." } # generates a string with @nbit dots: "..."
        dots = "(" + dots + ")"    # adds the parentesys: "(...)"
        rexp = Regexp.new( dots )  # converst dots into a regulare expression: /(...)/
        str_a = str.split( rexp )  # splits the string: eg."000111110011" -> ["","000","","111","","110","","011"]
        str_a = str_a.delete_if{ |v| v == ""}  # ["000", "111", "110", "011"]
        str_a.size.times do |i| 
            if i <= str_a.size-2
                # the first bit of the element with an odd index is the sign of the floating number
                str_a[ 2*i ][0..0] == "0" ? sign = "+" : sign = "-"
                # this returns: [ +00.7 , +2.6 ], each element in the array is in decimal codification
                cr << (sign+str_a[ 2*i ][1..-1].to_i(2).to_s+"."+str_a[2*i+1].to_i(2).to_s).to_f
            end # if
            return cr if cr.size == str_a.size / 2
        end # str_a.size.times
    end # decode
  end # class Optimizer
end # module GA#


if __FILE__ == $0
  # Test function
  f = lambda {|p| p[0]**2 + p[1]**2 } # a trivial parabola
  #f = lambda { |p| p[0] ** 2 - 4 * p[0] + p[1] ** 2  -p[1] - p[0] * p[1]}
  #f = lambda { |p| ( 1 - p[0] ) ** 2 + 100 * ( p[1] - p[0] ) ** 2 } # Rosenbroke function
  
  # Instantiate the optimizer, with tolerance and dimension
  opt = GA::Optimizer.new( :tol => 1,
      :p_mutation  => 0.2,
      :p_crossover => 0.8,
      :i_o         => { :X =>[5,10] , :Y=>[-10.23,5.234] },
      :npop        => 50,
      :ncr         => 200
    )
  opt.loop {|p| f.call(p)}
end
