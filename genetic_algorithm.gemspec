# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = %q{genetic_algorithm}
  s.version = "0.0"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Carlos Maximiliano Giorgio Bort"]
  s.date = %q{2011-03-16}
  s.description = %q{Simple genetic algorithm for functions minimization.}
  s.email = %q{maximiliano_giorgio@yahoo.it}
  s.files = %w[README.txt lib/genetic_algorithm.rb test.rb]
  s.has_rdoc = true
  s.homepage = %q{https://github.com/maximiliano1985}
  s.rdoc_options = ["--inline-source", "--charset=UTF-8"]
  s.require_paths = ["lib"]
  s.rubyforge_project = %q{genetic_algorithm}
  s.rubygems_version = %q{1.6.2}
  s.has_rdoc = true
  s.summary = %q{Simple genetic algorithm for function minimization. Needs gnuplot and the gem gnuplotr to work properly.}

end