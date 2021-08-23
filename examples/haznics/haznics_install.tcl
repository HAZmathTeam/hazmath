#!/usr/bin/wish

########################################################################
# HAZmath+FEniCS Python Interface (HAZniCS) 
#
# Tcl install script. Installs the HAZmath library, the swig library,
# cbc.block and FEniCS_ii and dependencies.
# 
########################################################################
proc InstallError { x } {
    puts stderr {*****************************************************************}
    puts stderr "ERROR: $x"
    puts stderr {*****************************************************************}
    exit
}
proc InstallStatus { s } {
    global status
    set status $s
    update idletasks
}

proc InstallHelper { s } {
    global helper
    global helper_add
    set h $s
    if { [llength $h] > 0 } {
	set helper [lindex $h 0]
    } else {
	set helper {}
    }
    if { [llength $h] > 1 } {
	set helper_add [lindex $h 1]
    } else {
	set helper_add {}
    }
    update idletasks
}

proc SetDirs { } {
    global dirs
    global dirlist
    global add_dirs
    global add_dirlist
    set haz_list [lreplace [file split [pwd]] end-1 end]
    set dirs(hazmath_dir) [file join {*}$haz_list]
    set dirs(haznics_dir) [pwd]
    set haznics_list [file split [pwd]]
    set add_dirs(haz_build_dir) [file join {*}$haz_list BUILD_HAZ]
    set add_dirs(quadpy_dir) [file join {*}$haznics_list quadpy]
    set add_dirs(fenics_ii_dir) [file join {*}$haznics_list fenics_ii]
    set add_dirs(cbcblock_dir) [file join {*}$haznics_list [join [list cbc . block] ""]]
    set add_dirs(ulfy_dir) [file join {*}$haznics_list ulfy]
    #[list hazmath_dir haznics_dir]
    set dirlist [array names dirs]
    #[list haz_build_dir quadpy_dir fenics_ii_dir cbcblock_dir ulfy_dir]
    set add_dirlist [array names add_dirs]
    update idletasks
}

proc Help {} {
    global prompthelp
    if {[info exists prompthelp(ok)]} {unset prompthelp(ok)}
    set f [toplevel .prompthelp -borderwidth 10 -width 30]
    set b [frame $f.buttons -bd 10]
    ScrolledText $f 96 24
    button $b.close -text close -command {set prompthelp(ok) 1}
    $f.text configure \
	-font -*-Helvetica-bold-r-*--*-120-* \
	-background white -foreground black -state normal
    label $f.label -relief flat\
	-text "HAZniCS install help" -border 0 \
	-font  -*-Helvetica-bold-r-*--*-120-* 
    pack $b.close -side left
    pack $f.label -side top
    pack $f.buttons -side top
    pack $f.yscroll -side right -fill y
    pack $f.label  -side top -fill x -anchor e
    pack $f.text  -side top -fill both -expand true
    pack $f.xscroll -side bottom -fill x
    pack $f.buttons -side bottom
    if {[catch {open "HAZniCS_GUI_help.txt" r} in1]} {
	set line "Could not read help file HAZniCS_GUI_help.txt"
	InstallHelper [list $line]
	$f.text insert end "$line\n"
    } else {
	while { [gets $in1 line] >=0 } {
	    #	    puts $line
	    $f.text insert end "$line\n"
	}
	close $in1
    }
    $f.text configure -state disabled
    grab $f
    tkwait variable prompthelp(ok)
    grab release $f
    if {$prompthelp(ok)} {destroy $f}
}

proc ScrolledText { parent width height } {
    # This proc follows closely the Brent's book
    catch {frame $parent}
    text $parent.text
    $parent.text configure -setgrid true \
	-xscrollcommand [list $parent.xscroll set] \
	-yscrollcommand [list $parent.yscroll set] \
	-wrap word -width $width -height $height
    scrollbar $parent.xscroll -orient horizontal \
	-command [list $parent.text xview]
    scrollbar $parent.yscroll -orient vertical \
	-command [list $parent.text yview]
}

proc ScrolledText1 { parent width height } {
    # This proc follows closely the Brent's book
    catch {frame $parent}
    text $parent.text 
    $parent.text configure -setgrid true \
	-xscrollcommand [list $parent.xscroll set] \
	-yscrollcommand [list $parent.yscroll set] \
	-wrap word -width $width -height $height \
	-borderwidth 5 -relief sunken
    scrollbar $parent.xscroll -orient horizontal \
	-command [list $parent.text xview]
    scrollbar $parent.yscroll -orient vertical \
	-command [list $parent.text yview]
    pack $parent.yscroll -side right -fill y
    pack $parent.xscroll -side bottom -fill x
}

proc ItemsSelect { w y } {
    #    $w select set anchor [$w nearest $y]
    set ix [$w nearest $y]
    $w select set $ix
    $w see $ix
}

proc  Uninstall { f } {
    global status
    global haz_build_dir
    global dirs
    global dirlist
    global add_dirs
    global add_dirlist
    SetDirs
    $f configure -state normal
    $f insert end "\nUninstall: Removing HAZmath install: lib and swig_files\n"
    update idletasks
    if {[file isdirectory $add_dirs(haz_build_dir)]} {
	set command0 "|make VERBOSE=1 -C $dirs(hazmath_dir) distclean"
	catch "open [list $command0] r" in
	update idletasks
	while {[gets $in line] >=0}  {
	    $f insert end "Uninstall: $line\n"
	    update idletasks
	    $f see end
	}
	close $in
    }  
    if {[info exists add_dirs(haz_build_dir)]} {unset add_dirs(haz_build_dir)}
    InstallHelper [list "Uninstall: HAZmath uninstalled"]
    update idletasks
    foreach el $add_dirlist {
	if {![info exists add_dirs($el)]} {
	    continue
	}
	if {[file isdirectory "$add_dirs($el)"]} {
	    if { [ catch [list exec /bin/rm -rf $add_dirs($el)]  msg ] } {
		$f insert end "\nERROR: COULD NOT REMOVE $add_dirs($el); $msg"	
	    } else {
		$f insert end "\nUninstall: Removed $add_dirs($el)"
	    }
	} else {
	    $f insert end "\nDirectory $add_dirs($el) does not exist or is not a directory!"
	}
	if {[info exists add_dirs($el)]} {unset add_dirs($el)}
	update idletasks
    }
    $f insert end "\n\n"
    $f see end
    $f configure -state disabled
    InstallHelper [list "" "Uninstall DONE"]
    InitOpts
}

proc Configure { f } {
    global build_type
    global pack_opts
    global packlist
    global opts
    global optlist
    global dirs
    global dirlist
    global add_dirs
    global add_dirlist
    global option_lbl
    global pack_lbl
    global status
    SetDirs
    InstallHelper [list "Configuring..."]
    #
    set make_flags [list "haznics=yes"]
    ##
    foreach el $optlist {
	if {$opts($el)} {
	    set  make_flags [lappend make_flags \
				 "$el=yes" ]
	} else {
	    #do nothing
	}
    }
    set xcommand(make)  "|/bin/sh -c  \"make -C $dirs(hazmath_dir) config $make_flags \" "
    set in  [eval {open} [list $xcommand(make) r]]
    $f configure -state normal
    while {[gets $in line] >=0}  {
	$f insert end "$line\n"
	update idletasks
	$f see end
    }
    $f see end
    $f configure -state disabled 
    InstallHelper [list "Configure...DONE"]
}
proc RunIt { f command0 } {
    global status
    global dirs
    global dirlist    
    global add_dirs
    global add_dirlist    
    $f configure -state normal
    if {![info exists add_dirs(haz_build_dir)]} {
    	InstallHelper [list "Run Configure First!"]
    	$f insert end "\nERROR: Build Dir NOT set. Run Configure First!\n"
    	update idletasks
    	InitOpts
    } elseif {![file isdirectory $add_dirs(haz_build_dir)]} {
    	InstallHelper [list "Build Dir does not exist; Run Configure First!"]
    	$f insert end "\nERROR: Build Dir does not exist or is not a directory; Run Configure First!\n"
    	InitOpts
    	update idletasks
    } else {	
	set flag 1
	if {![file exists $add_dirs(haz_build_dir)/Makefile]} {
	    InstallHelper [list "Run Configure first!"]
	    $f insert end "\nFile $add_dirs(haz_build_dir)/Makefile does not exist. Run Configure first!\n"
	    update idletasks
	} else {
	    InstallHelper [list "Installing... Please wait..."]
	    switch -exact -- [CheckInpData $command0] {
#		docs { set xcommand(make) "|make -C $dirs(hazmath_dir) $command0 2>@ stdout"}		
#		headers { set xcommand(make) "|make -C $dirs(hazmath_dir) $command0 2>@ stdout"}
		install { set xcommand(make) "|make -C $dirs(hazmath_dir) $command0 2>@ stdout"}
		default { set flag 0 } 
	    }
	    update idletasks
	    if { $flag == 1 } {
		catch "open [list $xcommand(make)] r" in
		# $f configure -state normal
		update idletasks
		while {[gets $in line] >=0}  {
		    $f insert end "$line\n"
		    update idletasks
		    $f see end
		}
		close $in
	    }
	    InstallHelper [list "Installed {HAZmath SWIG_libs}"]
	    update idletasks
	}
    }
    $f see end
    $f configure -state disabled 
}

proc InstallAdd { f } {
    global status
    global dirs
    global dirlist    
    global add_dirs
    global add_dirlist    
    global pack_opts
    global pack_repos
    global packlist
    $f configure -state normal
    set pinstalled {}
    foreach p $packlist {
	if { ! $pack_opts($p) } {
	    $f insert end "\nPackage $p is not to be installed\n"
	    continue
	} else {
	    set pinstalled [lappend pinstalled $p]
	}
	set key [join [list $p "_dir"] ""]	
	InstallHelper [list "" "Installing application packages... Please wait"]
	if {![file isdirectory $add_dirs($key)]} {
	    $f insert end "\nClonning $pack_repos($p) in $add_dirs($key)\n"
	    InstallHelper [list "" "Clonning $p"]
	    set command0 "|git clone $pack_repos($p) 2>@ stdout"
	    catch "open [list $command0] r" in
	    # $f configure -state normal
	    update idletasks
	    while {[gets $in line] >=0}  {
		$f insert end "$line\n"
		update idletasks
		$f see end
	    }
	    catch { close $in }
	}
	update idletasks
	set xdir [pwd]
	catch { cd $add_dirs($key) }
	if { $p == "quadpy" } {
	    catch { exec git checkout v0.12.10 }
	}
	set command0 "|python3 setup.py install --user  2>@ stdout"
	##    set command0 "|python3 setup.py install --user"
	set in1 [open $command0 r]
	while {[gets $in1 line] >=0}  {
	    $f insert end "$line\n"
	    update idletasks
	    $f see end
	}
	catch { cd $xdir }
	if { [catch { close $in1 } msg ]} {
	    InstallHelper [list "" "ERROR: $p could not be installed"]	
	    $f insert end "\nERROR in executing $command0\n"
	} else {
	    InstallHelper [list "" "$p is now installed"]
	    $f insert end "\n$p is now installed\n"
	}
    }
    $f insert end [concat "\nInstalled:" $pinstalled "\n"]
    InstallHelper [list "" [list "Installed:" $pinstalled]]
    $f see end
    $f configure -state disabled 
    update idletasks
}

proc GentleExit {} {
    exit
}

proc CheckInpData { x } {
    set y {}
    regsub "^\[ \t\]*" $x {} y
    regsub "\[ \t\]*\$" $y {} y
    regsub -all -- "\[ \t\]\[ \t\]*" $y { } y
    return $y
}
proc InitOpts { } {
    #Initialize options
    global pack_opts
    global packlist
    global pack_lbl
    global pack_repos
    global opts
    global optlist
    global option_lbl
    #
    global dirs
    global dirlist
    #
    global add_dirs
    global add_dirlist
    global status
    global build_type
    set option_lbl(verbose) "Verbose output \?"
    set option_lbl(shared) "Build shared library \?"
    set option_lbl(lapack) "Build with LAPACK support \?"
    set option_lbl(suitesparse) "Build with UMFPACK support \?"
    set option_lbl(debug) "Debugging \?"
    set build_type "RELEASE"    
    set optlist [list verbose shared suitesparse lapack debug]
    set opts(maxlen) -1
    foreach el $optlist {
	set opts($el) 1
	set len [expr [string length $option_lbl($el)] + 2]
	if { $opts(maxlen) < $len } {
	    set opts(maxlen) $len
	}
    }
    set opts(verbose) 0
    set opts(debug) 0
    ###################################################
    set pack_lbl(cbcblock) "Package CBC.BLOCK \?"
    set pack_lbl(fenics_ii) "Package FEniCS_ii \?"
    set pack_lbl(ulfy) "Package ULFY \?"
    set pack_lbl(quadpy) "Package QUADPY \?"
    set packlist [list cbcblock fenics_ii ulfy quadpy]
    # repos for the different packages:
    set pack_repos(cbcblock) "https://bitbucket.org/fenics-apps/cbc.block.git"
    set pack_repos(fenics_ii) "https://github.com/MiroK/fenics_ii.git"
    set pack_repos(ulfy) "https://github.com/MiroK/ulfy.git"
    set pack_repos(quadpy) "https://github.com/nschloe/quadpy.git"
    ###################################
    set packlist [list cbcblock quadpy fenics_ii ulfy]
    set pack_opts(maxlen) -1
    foreach el $packlist {
	set pack_opts($el) 1
	set len [expr [string length $pack_lbl($el)] + 2]
	if { $pack_opts(maxlen) < $len } {
	    set pack_opts(maxlen) $len
	}
    }
    SetDirs
    InstallHelper [list "Initialization is DONE..."]
}

proc ShowOpts { f flag } {
    ## if flag true then we do hazmath otherwise we do the opts
    global pack_opts
    global packlist
    global pack_lbl
    global opts
    global optlist
    global option_lbl
    #
    global dirs
    global dirlist
    #
    global add_dirs
    global add_dirlist
    global compiler_flags
    global status
    global build_type
    if { $flag } {
	set maxlen $opts(maxlen)
	unset opts(maxlen)
	foreach el $optlist {
	    frame $f.$el
	    label $f.$el.l -relief flat -width $maxlen -anchor w \
		-text $option_lbl($el) -border 0 \
		-font  -*-Helvetica-bold-r-*--*-140-* 
	    radiobutton $f.$el.y -variable opts($el) \
		-text "Yes" -value 1 -state normal
	    radiobutton $f.$el.n -variable opts($el) \
		-text "No" -value 0 -state normal 
	    pack $f.$el.l -side left 
	    pack $f.$el.n -side right -fill both
	    pack $f.$el.y -side right -fill both
	    pack $f.$el -side top -fill both
	}
	set opts(maxlen) $maxlen
	update idletasks
    } else {
	set maxlen $pack_opts(maxlen)
	unset pack_opts(maxlen)
	foreach el $packlist {
##	    puts [concat $el $pack_opts($el)]
	    frame $f.$el
	    label $f.$el.l -relief flat -width $maxlen -anchor w \
		-text $pack_lbl($el) -border 0 \
		-font  -*-Helvetica-bold-r-*--*-140-* 
	    radiobutton $f.$el.y -variable pack_opts($el) \
		-text "Yes" -value 1 -state normal
	    radiobutton $f.$el.n -variable pack_opts($el) \
		-text "No" -value 0 -state normal 
	    pack $f.$el.l -side left 
	    pack $f.$el.n -side right -fill both
	    pack $f.$el.y -side right -fill both
	    pack $f.$el -side top -fill both
	}	
	set pack_opts(maxlen) $maxlen
	update idletasks
    }
}
#########################################################
global opts
global optlist
global dirs
global dirlist
global add_dirs
global add_dirlist
global status
global helper
global helper_add

InitOpts 

frame .f

frame .f.c -borderwidth 10
frame .f.c.haz -relief groove -borderwidth 10 
frame .f.c.add -relief groove -borderwidth 10 
frame .f.c.haz.buttons -borderwidth 10
frame .f.c.add.buttons -borderwidth 10

frame .f.images -background "#00008B"
pack .f.images -side top -fill x
label .f.images.logo -background "#00008B" -text "HAZniCS" \
    -background "#00008B" -foreground yellow -font { Helvetica -24 bold } 
pack .f.images.logo -side left
label .f.images.text -text  "HAZmath+FEniCS Python Interface (HAZniCS)\nInstallation script" \
    -background "#00008B" -foreground "#FFFFFF" -font { Helvetica -18 bold } 
pack .f.images.text 

pack .f  -expand true -fill both

pack .f.c -side left -fill both
pack .f.c.add -side top -fill x 
pack .f.c.add.buttons -side bottom -fill x 
pack .f.c.haz -side bottom  -fill x 
pack .f.c.haz.buttons -side bottom -fill x 

##########main 

ScrolledText1 .f 96 24

## quit, uninstall help buttons. 
frame .f.buttons -borderwidth 10 
button .f.buttons.quit -text "Quit" -font { Helvetica -14 bold }  -background "#00008b" -foreground "#FFFFFF" -command {GentleExit}
button .f.buttons.help -text Help -font { Helvetica -14 bold }  -background "#008b00" -foreground "#FFFFFF" -command {Help}
button .f.buttons.uninstall -text "Uninstall all"  -font { Helvetica -14 bold }  -background "#8b0000" -foreground "#FFFFFF"  -command {Uninstall .f.text}
## command window. 
pack .f.buttons.help -side right
pack .f.buttons.quit -side right
pack .f.buttons.uninstall -side left
pack .f.buttons -side top  -fill x
update idletasks

## Text window:
pack .f.text -side right -fill both
.f.text configure -state disabled

#   grab release .f.text
button .f.c.haz.buttons.config -text "Configure HAZniCS" -command {Configure .f.text}
button .f.c.haz.buttons.install -text "Install HAZniCS" -command {RunIt .f.text install}

pack .f.c.haz.buttons.config .f.c.haz.buttons.install  -side left -fill x

button .f.c.add.buttons.install -text "Install Packages" -command {InstallAdd .f.text }

pack .f.c.add.buttons.install  -side left -fill x

frame .f.c.add.opts -borderwidth 10

frame .f.c.haz.opts -borderwidth 10

label .f.c.haz.opts.title -textvar title -relief raised \
    -background "#000000" -foreground "#FFFFFF" -font { Helvetica -14 bold } 

label .f.c.add.opts.title -textvar title_add -relief raised \
    -background "#000000" -foreground "#FFFFFF" -font { Helvetica -14 bold } 

set title "Install HAZniCS libraries"
set title_add "Install Application Packages"

##label .f.opts.lstatus -textvar status -relief raised -width 35  

pack .f.c.haz.opts.title -side top -fill both
pack .f.c.add.opts.title -side top -fill both

label .f.c.haz.opts.lhelper -textvar helper -relief flat -width 55 \
    -background "#00008B" -foreground "yellow" -font { Helvetica -14 bold } 

label .f.c.add.opts.lhelper -textvar helper_add -relief flat -width 55 \
    -background "#00008B" -foreground "#FFFFFF" -font { Helvetica -14 bold } 

##label .f.opts.lstatus -textvar status -relief raised -width 35  

pack .f.c.haz.opts.lhelper -side bottom -fill both
pack .f.c.add.opts.lhelper -side bottom -fill both
##pack .f.opts.lstatus -side bottom 



ShowOpts .f.c.haz.opts true

#EntryList .f.c.haz.opts

pack .f.c.haz.opts -side top -fill both

ShowOpts .f.c.add.opts false

#EntryList .f.c.add.opts

pack .f.c.add.opts -side bottom -fill both




