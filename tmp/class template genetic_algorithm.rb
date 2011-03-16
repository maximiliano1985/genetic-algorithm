#!/usr/bin/env ruby
#
#  Created by Carlos Maximiliano Giorgio Bort on 2011-03-06.
#  Copyright (c) 2011 University of Trento. All rights reserved.
#

require 'rubygems'
require 'matrix'
require 'gnuplotr'

class Vector
  if RUBY_VERSION.split(".").join.to_i < 190
    # Add method for division under Ruby 1.8.x
    def /(other); self * (1.0 / other); end
  end
  # More compact inspect version
  def inspect; "V[#{self.to_a * ","}]"; end
end

module GA

  class Population
    attr_reader :dimension, :size

    def initialize(dim, size)
        @dimension = dim # is number of cromosomes in the population
        @size = size     # is the problem size
    end
    
    def encode
    end
    
    def decode
    end
    
    def performance
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
        :dim   => 2,
        :exp_f => 1.5, 
        :cnt_f => 0.5,
        :tol   => 0.001,
        :niter => 50,
        :pconv => true,
        :plotopt=> {:title  => 'Nelder Mean Method Convergence',
                     :xlabel => 'No. iteration',
                     :ylabel => 'Objective function value',
                     :yrange => [ -1 , 3 ],
                     :grid   => "set grid"
                     }
      }
      @cfg.merge! args
      @simplex = Simplex.new(@cfg[:dim])
      @start_points = []
      @status = :filling
      @iteration = 0
      if @cfg[:pconv] == true # this is the plot
          @gp = GNUPlotr.new("/opt/local/bin/gnuplot")
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
      end
    end
     
    def evolve
    end
    
    def mutation
    end
    
    def crossover
    end
    
    def invertion
    end
    
    
    # Sets the starting points of the +Simplex+.
    # @param [Array] ary an +Array+ of +Vector+s
    # @raise [ArgumentError] unless +ary+ is an +Array+ of +@cfg[:dim]+ +Vector+s
    def start_points=(ary)
      raise ArgumentError, "As start_points I need an Array instead of a #{ary.class}" unless ary.kind_of? Array
      raise ArgumentError, "Array size must be #{@cfg[:dim]} but it's #{ary.size}" unless ary.size == @cfg[:dim]
      content = true
      ary.each { |v|  content = false unless v.kind_of? Vector}
      raise ArgumentError, "Array elements must be Vectors" unless content
      @start_points = ary
    end
    
    # Adds a new point to the +Simplex+.
    # @param [Array <Vector,Numeric>] point a new point with its value.
    # @raise [ArgumentError] is point is not of the expected type
    # @return [Array <Vector,Numeric>] the new point
    def <<(point)
      raise ArgumentError, "Need an Array" unless point.size == 2
      raise ArgumentError, "point[0] must be Vector" unless point[0].kind_of? Vector 
      raise ArgumentError, "point[1] must be Numeric" unless point[1].kind_of? Numeric
      @simplex[point[0]] = point[1]
    end
    
    # Starts the minimization loop. It expects the block containing the 
    # function evaluation.
    # @example
    #   f = lambda {|point| point[0]**2 + point[1]**2}
    #   opt = NMM::Optimizer.new( :dim => 3, :tol => 1E-5)
    #   opt.start_points = [Vector[10,37],Vector[7,2],Vector[51,32]]
    #   opt.loop {|point| f.call(p)}
    # @param [Array] ary if given, it gets filled with the optimization path
    # @yield [point] the block must return the evaluation of the function in +point+
    # @yieldparam [Vector] point the point which the function has to be evaluated at
    # @yieldreturn [Float] the evaluation
    # @raise [ArgumentError] unless a block is given
    def loop(ary = nil)
      raise ArgumentError, "Block needed" unless block_given?
      fx = nil
      until converged? or @iteration > @cfg[:niter] do 
        next_point = step(fx)
        unless next_point[1] then # if f(next_point) is "nil" do this...
          puts "Reflecting at #{next_point}"
          fx = yield(next_point[0])
        else # if f(next_point) is a number (or better, it isn't "nil"), do...
          next_point[1] = yield(next_point[0]) if next_point[1] == 0
          self << next_point
          fx = nil
          if @status != :filling
            puts "Running the #{@iteration}th iteration" unless @iteration == @cfg[:niter]
            puts "Reached max number of iteration [#{@iteration}]" if @iteration == @cfg[:niter]
            
            # these lines ar used to do a convergence plot
            if @cfg[:pconv] == true
                if @iteration == 0 # at first iteration initialize the matrix containing simplex data
                    it = []
                    sdata = []
                end
                it << @iteration
                sdata << [self.simplex[:l][1], self.simplex[:g][1], self.simplex[:h][1], self.simplex.norm ]
                sleep 0.1
                # a. initialize the data sets for the plot
                @gp.new_series(:simplex_l)
                @gp.new_series(:simplex_g)
                @gp.new_series(:simplex_h)
                @gp.new_series(:simplex_norm)
                # b. fill the data sets
                it.each do |i| 
                   @gp.series[:simplex_l]    << [ i, sdata[i][0] ]
                   @gp.series[:simplex_g]    << [ i, sdata[i][1] ]
                   @gp.series[:simplex_h]    << [ i, sdata[i][2] ]
                   @gp.series[:simplex_norm] << [ i, sdata[i][3] ] 
                end
                # c. close the data sets
                @gp.series[:simplex_l].close
                @gp.series[:simplex_g].close
                @gp.series[:simplex_h].close                
                @gp.series[:simplex_norm].close
                # d. plot the data sets
                @gp.plot :simplex_l      , "with lines"
                @gp.replot :simplex_h    , "with lines"
                @gp.replot :simplex_g    , "with lines"
                @gp.replot :simplex_norm , "with lines linewidth 3"
            end # if @cfg
            p self.status
            log
            @iteration += 1
          end # if @status  
        end # unless next_point
        ary << [next_point, @status] if ary.kind_of? Array    
        
        # when the solution converges, stores the plot data
        puts @gp.dump_input if @cfg[:pconv] == true && (converged? || @iteration > @cfg[:niter])
      end # until converged?
    end # def loop
    
    # This method prints out the track of the optimization.
    # @note You should override this method if you want to have a different 
    #   printout of minimization path.
    def log
      values = [ @simplex.points[-1][1],(@simplex.norm || 0)].flatten
      puts "New point at:\n #{@simplex.points[-1][0].to_a.inspect} -> %9.5f ||%9.5f||" % values 
      puts "______________________________________"
    end
    
    # Returns if the loop has converged, i.e. if the norm is less than +@cfg[:tol]+
    # @return [Bool] true if the loop has converged
    def converged?
      n = @simplex.norm
      if n then
        @simplex.norm < @cfg[:tol] # this returns true or false
      else
        false
      end
    end
    
    private
    def step(vr = nil)
      # Filling starting Simplex
      if @start_points.size > 0 then
        @status = :filling
        return [@start_points.shift, 0] 
      end
      
      # Reflecting Simplex
      unless vr
        @status = :reflecting
        return [@simplex[:r], nil]
      end
      x_new = [@simplex[:r], vr]
  
      # Checking Expansion/Contraction
      case
        when vr < @simplex[:l][1] # if the new point is lower than the lowest point in the simplex...r
          # Expansion
          @status = :expansion
          x_new = [ @simplex[:r]*@cfg[:exp_f] - @simplex[:c], 0 ]
      when vr >= @simplex[:h][1]
          # Contraction
          @status = :contraction1
          x_new = [ (@simplex[:h][0] + @simplex[:c])*@cfg[:cnt_f] , 0 ]
      when @simplex[:g][1] < vr && vr < @simplex[:h][1]
          # Contraction
          @status = :contraction2
          x_new = [(@simplex[:h][0] + @simplex[:r])*@cfg[:cnt_f], 0]
      end
      return x_new
    end
  end # class Optimizer
end # module NMM




if __FILE__ == $0 then
  # Test function
  #f = lambda {|p| p[0]**2 + 3*p[1]**2+10} # a trivial parabola
  #f = lambda { |p| p[0] ** 2 - 4 * p[0] + p[1] ** 2  -p[1] - p[0] * p[1]} # a non trivial parabola
  f = lambda { |p| ( 1 - p[0] ) ** 2 + 100 * ( p[1] - p[0] ) ** 2 } # the non trivial Rosenbroke function
  # Instantiate the optimizer, with tolerance and dimension (it is the dimension
  # of the simplex, so number of parameters + 1 !!!)
  opt = NMM::Optimizer.new( :dim => 3, :tol => 1E-5,:niter => 50, :exp_f => 2,:cnt_f => 0.5)
  
  # Define the starting points, i.e. the first guess simplex
  opt.start_points = [
    Vector[0.1,0.1],
    Vector[1.2,0.0],
    Vector[0.0,0.8] 
  ]
  
  # Start the loop, passing a block that evaluates the function at the point
  # represented by the block parameter p
  # For different logging, just override the Optimizer#log method
  opt.loop {|p| f.call(p)}
  
  # Final simplex configuration is available at the end.
  #p opt.simplex
  p opt.simplex

  puts "__________________________________________________"

  puts "The residual is #{opt.simplex.norm}, the number of iterations made are #{opt.iteration}"
  puts "The minimum is at #{opt.simplex[:l][0]} and his value is #{opt.simplex[:l][1]}"
  puts "__________________________________________________"
  # In order to save track of the solution, pass an existing array to the 
  # Optimizer#loop method:
  # ary = []
  # opt.loop(ary) {|p| f.call(p)}
  # puts "***"
  # p ary
end
