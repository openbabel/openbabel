#!/usr/bin/perl

#
# bin2hex.pl by Chami.com
# http://www.chami.com/tips/
#

# number of characters per line
$chars_per_line = 15;

# -------------------------------------

# language id
#
# 0 = Perl (default)
# 1 = C / C++
# 2 = Pascal / Delphi
#

$rem_begin  = "begin binary data:";
$rem_end    = "end binary data.";

# C / C++
$_var       = "\n\n/* $rem_begin */\n".
  "static const char $ARGV[1]"."[] = ".
  "/* %d */\n";
$_begin     = "{";
$_end       = "};\n";
$_break     = "\n";
$_format    = "0x%02X";
$_separator = ",";
$_comment   = "/* $rem_end ".
  "size = %d bytes */\n";

if(open(F, "awk 'substr(\$1,0,1) != \"#\" && NF > 0' ".$ARGV[0]." |"))
  {
    binmode(F);
    
#   $s = '';
    $i = 0;
    $count = 0;
    $first = 1;
    $s .= $_begin;
    while(!eof(F))
      {
	if($i >= $chars_per_line)
	  {
	    $s .= $_break;
	    $i = 0;
	  }
	if(!$first)
	  {
	    $s .= $_separator;
	  }
	$s .= sprintf($_format, ord(getc(F)));
	++$i;
	++$count;
	$first = 0;
      }

    if($i >= $chars_per_line)
    {
	$s .= $_break;
	$i = 0;
    }
    if(!$first)
    {
	$s .= $_separator;
    }
    $s .= "0x00";
    ++$count;

    $s .= "\n".$_end;
    $s .= sprintf $_comment, $count;
    $s .= "\n\n";
    
    $s = "\n".sprintf($_var, $count).$s;
    
    print $s;
    
    close( F );
  }
else
  {
    print
      "bin2hex.pl by Chami.com\n".
	"\n".
	  "usage:\n".
	    "  perl bin2hex.pl ".
	      " \n".
		"\n".
		  "   : path to the ".
		    "binary file\n".
		      "   : 0 = Perl, ".
			"1 = C/C++/Java, ".
			  "2 = Pascal/Delphi\n".
			    "\n";
  }

