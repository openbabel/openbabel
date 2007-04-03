# Creates a makefile to build the Open Babel-Ruby extension.

# This tells mkmf where to look for the Open Babel headers. There's probably a better way to do it, but the mkmf documentation doesn't say how.
if (ARGV.size == 0)
  ARGV << "--with-openbabel-include=../../include"
end

require 'mkmf'

dir_config('openbabel')

if have_library('openbabel')
  create_makefile('openbabel')
else
  puts "Install Open Babel first. If you've already compiled and installed Open Babel, you may need to run ldconfig."
end


