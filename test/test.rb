#!/usr/bin/env ruby
#  ga_my
#
#  Created by Carlos Maximiliano Giorgio Bort on 2011-03-06.
#  Copyright (c) 2011 University of Trento. All rights reserved.
#

require 'rubygems'
require 'gnuplotr'
require '../lib/genetic_algorithm'
# Test function
f = lambda {|p| p[0]**2 + p[1]**2 } # a trivial parabola
  
# Instantiate the optimizer, with tolerance and dimension
opt = GA::Optimizer.new( :tol => 1E-3,
    :p_mutation  => 0.2,
    :p_crossover => 0.8,
    :i_o         => { :X =>[-5,10] , :Y=>[-10.23,5.234] },
    :npop        => 50,
    :ncr         => 150
)
opt.loop {|p| f.call(p)}
