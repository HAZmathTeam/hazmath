# generate prototypes for Samba C code
# originally writen by A. Tridge, June 1996
# this modified version was taken from rsync 2.4.6

# removed some things that I do not need in SiPSMG (ltz, 2009)
# modified by JAdler for FEM_MHD (2010-May-07)

BEGIN {
  inheader=0;
  print "/*******************************************************************/  "
  print "/* This header file was automatically generated with \"make fheaders\".   */" 
  print "/* WARNING: DO NOT EDIT!!!                               */  "
  print "/*******************************************************************/  "
  print ""
  print "// Standard Includes"
  print "#include <stdlib.h>"
  print "#include <stdio.h>"
  print "#include <string.h>"
  print "#include <math.h>"
  print "#include <limits.h>"
  print "#include <getopt.h>"
  print "#include <sys/types.h>"
  print "#include <time.h>"
  print "#include <unistd.h>"
  print "// Internal Includes"
  print "#include \"macro.h\""
  print "#include \"grid.h\""
  print "#include \"sparse.h\""
  print "#include \"vec.h\""
  print "#include \"fem.h\""
  print "#include \"solver.h\""
  print "#include \"param.h\""
}

{
  if (inheader) {
    if (match($0,"[)][ \t]*$") || match($0,"[)][ \t][{][ \t]*$")) {
      inheader = 0;
      printf "%s;\n",$0;
    } else {
      printf "%s\n",$0;
    }
    next;
  }
}

/\/*! \\file/ {
    printf "\n/* In file: %s */\n",$3;  
}

/^static|^extern/ || !/^[a-zA-Z]/ || /[;]/ {
  next;
}

!/^INT|^REAL|^coordinates|^qcoordinates|^FILE|^OFF_T|^size_t|^off_t|^pid_t|^unsigned|^mode_t|^DIR|^user|^int|^char|^uint|^struct|^BOOL|^void|^double|^time|^dCSRmat|^dvector|^iCSRmat|^ivector|^AMG_data/ {
  next;
}

/[(].*[)][ \t]*$/ {
    printf "%s;\n",$0;
    next;
}

/[(]/ {
  inheader=1;
  printf "%s\n",$0;
  next;
}

END {
  print "\n/* End of header file */"
}
