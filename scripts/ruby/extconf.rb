# Creates a makefile to build the Open Babel-Ruby extension.

# Compensate for the fact that Ruby will try to build universal
# binaries on OS X by default
require 'rbconfig'
if Config::CONFIG["arch"] =~ /universal-darwin/
  ENV['ARCHFLAGS'] = case `uname -smr`.chomp
    when "i386" then '-arch i386'
    when "ppc"  then '-arch ppc'
  end
end

require 'mkmf'

dir_config('openbabel')

# Find a trivial header in order to add the proper include path
# to the build flags.
here = File.dirname(__FILE__)
find_header('inchi_api.h', '/usr/include/inchi', '/usr/include', here + '/../../include')

# Prevent Ruby 1.8.x from trying to compile and link the extension
# using gcc.
if RUBY_VERSION < "1.9"
  cxx = ''
  begin
    File.open('../Makefile', 'r').each_line do |line|
      if line =~ /CXX = /
        cxx = Regexp.last_match.post_match.chomp
      end
    end
  rescue Errno::ENOENT
    puts 'Please configure Open Babel before compiling the Ruby extension'
  end
  cpp_command(cxx)
end

if have_library('openbabel')
  with_ldflags("#$LDFLAGS -dynamic -flat_namespace") do #Enables cc to handle linking better.
  create_makefile('openbabel')
end
else
  puts "Install Open Babel first. If you've already compiled and installed Open Babel, you may need to run ldconfig."
end
