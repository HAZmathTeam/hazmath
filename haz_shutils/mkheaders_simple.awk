# generate prototypes for Samba C code
# originally writen by A. Tridge, June 1996
# this modified version was taken from rsync 2.4.6

# removed some things that I do not need in SiPSMG (ltz, 2009)
# modified by JAdler for FEM_MHD (2010-May-07)

BEGIN {
  inheader=0;
  print "/*******************************************************************/  "
  print "/* This header file was automatically generated with \"make headers\".   */"
  print "/* WARNING: DO NOT EDIT!!!                               */  "
  print "/*******************************************************************/  "
  print ""
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

!/^INT|^REAL|^coordinates|^mesh_struct|^qcoordinates|^FILE|^OFF_T|^size_t|^off_t|^pid_t|^unsigned|^mode_t|^DIR|^user|^int|^char|^uint|^struct|^SHORT|^BOOL|^void|^double|^time|^dCSRmat|^dvector|^iCSRmat|^ivector|^dCOOmat|^dDENSEmat|^iDENSEmat|^block_dCSRmat|^AMG_data|^AMG_param|^scomplex|^MG_blk_data|^HX_curl_data|^HX_div_data|^precond_block_data|^precond_data|^precond_ra_data|^precond|^smoother_data|^smoother_matvec|^PyObject|^subscomplex|^macrocomplex|^unigrid|^cube2simp|^input_grid|^coordsystem|^features|^locdetails/ {

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

