#!/usr/bin/env ruby
#  binary_cod
#
#  Created by Carlos Maximiliano Giorgio Bort on 2011-03-10
#  Copyright (c) 2011 University of Trento. All rights reserved.
#

class Integer
    # size is the number of bits used 
    def to_bin(size)
        num = self.to_i # the integer part
        bin = num.to_s(2)
        bin.size <= size ? inc = (size - bin.size) : (raise "Use more bits for a proper binary convertion")
        inc.times { |i| bin = "0" +  bin }
        return bin
    end
end
class String    
    def to_dec
        self.to_i(2)
    end
end

num = 4.2
num_i = num.to_i # the integer part
num_d = ( (num-num.to_i)*1000 ).to_i# the decimal part
size = [ num_i.to_s(2).size , num_d.to_s(2).size ].max
bin_i = (num_i).to_bin(size)
bin_d = (num_d).to_bin(size)
puts "#{num_i} -> #{bin_i}"
puts "#{num_d} -> #{bin_d}"
puts "#{num_i}.#{num_d} -> #{bin_i}.#{bin_d}"
# la codifica va fatta dinamicamente, cambiando "size" ad ogni loop
n_i = bin_i.to_dec
n_d = bin_d.to_dec
puts "#{bin_i}.#{bin_d} -> #{n_i}.#{n_d}"

