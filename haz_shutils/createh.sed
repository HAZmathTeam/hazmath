#!/bin/sed -f

/^      subroutine.*)/{
s/)/);/g
t ends
}
/^      subroutine/{
:begin
N
s/)/);/g
t ends
b begin
}
/^      subroutine/!{
d
}
:ends
/^      subroutine/{
s/\n     [^[:space:]]//g
s/\n[cC].*\n//g
s/^[cC].*$//g
s/\!.*\n//g
s/\!.*$//g
s/\n//g
s/ *//g
s/,[abcdefghopqrstuvwxyzABCDEFGHOPQRSTUVWXYZ]/REAL *&/g
s/REAL \*,/,REAL \*/g
s/([abcdefghopqrstuvwxyzABCDEFGHOPQRSTUVWXYZ]/REAL *&/g
s/REAL \*(/(REAL \*/g
s/,[ijklmnIJKLMN]/INT *&/g
s/INT \*,/,INT \*/g
s/([ijklmnIJKLMN]/INT *&/g
s/INT \*(/(INT \*/g
s/(/_(/g
s/subroutine/extern void /g
}

