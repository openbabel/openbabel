require 'mkmf'

dir_config('openbabel')
dir_config('openbabel-2.0')
have_library('openbabel')

create_makefile('openbabel')
