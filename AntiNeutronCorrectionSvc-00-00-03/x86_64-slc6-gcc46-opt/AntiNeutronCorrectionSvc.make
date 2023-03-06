#-- start of make_header -----------------

#====================================
#  Library AntiNeutronCorrectionSvc
#
#   Generated Sat Feb 19 05:16:29 2022  by lliu
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_AntiNeutronCorrectionSvc_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_AntiNeutronCorrectionSvc_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_AntiNeutronCorrectionSvc

AntiNeutronCorrectionSvc_tag = $(tag)

#cmt_local_tagfile_AntiNeutronCorrectionSvc = $(AntiNeutronCorrectionSvc_tag)_AntiNeutronCorrectionSvc.make
cmt_local_tagfile_AntiNeutronCorrectionSvc = $(bin)$(AntiNeutronCorrectionSvc_tag)_AntiNeutronCorrectionSvc.make

else

tags      = $(tag),$(CMTEXTRATAGS)

AntiNeutronCorrectionSvc_tag = $(tag)

#cmt_local_tagfile_AntiNeutronCorrectionSvc = $(AntiNeutronCorrectionSvc_tag).make
cmt_local_tagfile_AntiNeutronCorrectionSvc = $(bin)$(AntiNeutronCorrectionSvc_tag).make

endif

include $(cmt_local_tagfile_AntiNeutronCorrectionSvc)
#-include $(cmt_local_tagfile_AntiNeutronCorrectionSvc)

ifdef cmt_AntiNeutronCorrectionSvc_has_target_tag

cmt_final_setup_AntiNeutronCorrectionSvc = $(bin)setup_AntiNeutronCorrectionSvc.make
cmt_dependencies_in_AntiNeutronCorrectionSvc = $(bin)dependencies_AntiNeutronCorrectionSvc.in
#cmt_final_setup_AntiNeutronCorrectionSvc = $(bin)AntiNeutronCorrectionSvc_AntiNeutronCorrectionSvcsetup.make
cmt_local_AntiNeutronCorrectionSvc_makefile = $(bin)AntiNeutronCorrectionSvc.make

else

cmt_final_setup_AntiNeutronCorrectionSvc = $(bin)setup.make
cmt_dependencies_in_AntiNeutronCorrectionSvc = $(bin)dependencies.in
#cmt_final_setup_AntiNeutronCorrectionSvc = $(bin)AntiNeutronCorrectionSvcsetup.make
cmt_local_AntiNeutronCorrectionSvc_makefile = $(bin)AntiNeutronCorrectionSvc.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)AntiNeutronCorrectionSvcsetup.make

#AntiNeutronCorrectionSvc :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'AntiNeutronCorrectionSvc'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = AntiNeutronCorrectionSvc/
#AntiNeutronCorrectionSvc::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

AntiNeutronCorrectionSvclibname   = $(bin)$(library_prefix)AntiNeutronCorrectionSvc$(library_suffix)
AntiNeutronCorrectionSvclib       = $(AntiNeutronCorrectionSvclibname).a
AntiNeutronCorrectionSvcstamp     = $(bin)AntiNeutronCorrectionSvc.stamp
AntiNeutronCorrectionSvcshstamp   = $(bin)AntiNeutronCorrectionSvc.shstamp

AntiNeutronCorrectionSvc :: dirs  AntiNeutronCorrectionSvcLIB
	$(echo) "AntiNeutronCorrectionSvc ok"

#-- end of libary_header ----------------

AntiNeutronCorrectionSvcLIB :: $(AntiNeutronCorrectionSvclib) $(AntiNeutronCorrectionSvcshstamp)
	@/bin/echo "------> AntiNeutronCorrectionSvc : library ok"

$(AntiNeutronCorrectionSvclib) :: $(bin)AN_XYZ_ErrorMatrx.o $(bin)AntiNeutronCorrectionSvc.o $(bin)AntiNeutronTrk.o $(bin)AntiNeutronCorrectionSvc_entries.o $(bin)AntiNeutronCorrectionSvc_load.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(AntiNeutronCorrectionSvclib) $?
	$(lib_silent) $(ranlib) $(AntiNeutronCorrectionSvclib)
	$(lib_silent) cat /dev/null >$(AntiNeutronCorrectionSvcstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(AntiNeutronCorrectionSvclibname).$(shlibsuffix) :: $(AntiNeutronCorrectionSvclib) $(AntiNeutronCorrectionSvcstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" AntiNeutronCorrectionSvc $(AntiNeutronCorrectionSvc_shlibflags)

$(AntiNeutronCorrectionSvcshstamp) :: $(AntiNeutronCorrectionSvclibname).$(shlibsuffix)
	@if test -f $(AntiNeutronCorrectionSvclibname).$(shlibsuffix) ; then cat /dev/null >$(AntiNeutronCorrectionSvcshstamp) ; fi

AntiNeutronCorrectionSvcclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)AN_XYZ_ErrorMatrx.o $(bin)AntiNeutronCorrectionSvc.o $(bin)AntiNeutronTrk.o $(bin)AntiNeutronCorrectionSvc_entries.o $(bin)AntiNeutronCorrectionSvc_load.o

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

ifeq ($(INSTALLAREA),)
installarea = $(CMTINSTALLAREA)
else
ifeq ($(findstring `,$(INSTALLAREA)),`)
installarea = $(shell $(subst `,, $(INSTALLAREA)))
else
installarea = $(INSTALLAREA)
endif
endif

install_dir = ${installarea}/${CMTCONFIG}/lib
AntiNeutronCorrectionSvcinstallname = $(library_prefix)AntiNeutronCorrectionSvc$(library_suffix).$(shlibsuffix)

AntiNeutronCorrectionSvc :: AntiNeutronCorrectionSvcinstall

install :: AntiNeutronCorrectionSvcinstall

AntiNeutronCorrectionSvcinstall :: $(install_dir)/$(AntiNeutronCorrectionSvcinstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(AntiNeutronCorrectionSvcinstallname) :: $(bin)$(AntiNeutronCorrectionSvcinstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(AntiNeutronCorrectionSvcinstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(AntiNeutronCorrectionSvcinstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(AntiNeutronCorrectionSvcinstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(AntiNeutronCorrectionSvcinstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(AntiNeutronCorrectionSvcinstallname) $(install_dir)/$(AntiNeutronCorrectionSvcinstallname); \
	      echo `pwd`/$(AntiNeutronCorrectionSvcinstallname) >$(install_dir)/$(AntiNeutronCorrectionSvcinstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(AntiNeutronCorrectionSvcinstallname), no installation directory specified"; \
	  fi; \
	fi

AntiNeutronCorrectionSvcclean :: AntiNeutronCorrectionSvcuninstall

uninstall :: AntiNeutronCorrectionSvcuninstall

AntiNeutronCorrectionSvcuninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(AntiNeutronCorrectionSvcinstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(AntiNeutronCorrectionSvcinstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(AntiNeutronCorrectionSvcinstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(AntiNeutronCorrectionSvcinstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),AntiNeutronCorrectionSvcclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)AN_XYZ_ErrorMatrx.d

$(bin)$(binobj)AN_XYZ_ErrorMatrx.d :

$(bin)$(binobj)AN_XYZ_ErrorMatrx.o : $(cmt_final_setup_AntiNeutronCorrectionSvc)

$(bin)$(binobj)AN_XYZ_ErrorMatrx.o : $(src)AN_XYZ_ErrorMatrx.cxx
	$(cpp_echo) $(src)AN_XYZ_ErrorMatrx.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(lib_AntiNeutronCorrectionSvc_pp_cppflags) $(AN_XYZ_ErrorMatrx_pp_cppflags) $(use_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(lib_AntiNeutronCorrectionSvc_cppflags) $(AN_XYZ_ErrorMatrx_cppflags) $(AN_XYZ_ErrorMatrx_cxx_cppflags)  $(src)AN_XYZ_ErrorMatrx.cxx
endif
endif

else
$(bin)AntiNeutronCorrectionSvc_dependencies.make : $(AN_XYZ_ErrorMatrx_cxx_dependencies)

$(bin)AntiNeutronCorrectionSvc_dependencies.make : $(src)AN_XYZ_ErrorMatrx.cxx

$(bin)$(binobj)AN_XYZ_ErrorMatrx.o : $(AN_XYZ_ErrorMatrx_cxx_dependencies)
	$(cpp_echo) $(src)AN_XYZ_ErrorMatrx.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(lib_AntiNeutronCorrectionSvc_pp_cppflags) $(AN_XYZ_ErrorMatrx_pp_cppflags) $(use_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(lib_AntiNeutronCorrectionSvc_cppflags) $(AN_XYZ_ErrorMatrx_cppflags) $(AN_XYZ_ErrorMatrx_cxx_cppflags)  $(src)AN_XYZ_ErrorMatrx.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),AntiNeutronCorrectionSvcclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)AntiNeutronCorrectionSvc.d

$(bin)$(binobj)AntiNeutronCorrectionSvc.d :

$(bin)$(binobj)AntiNeutronCorrectionSvc.o : $(cmt_final_setup_AntiNeutronCorrectionSvc)

$(bin)$(binobj)AntiNeutronCorrectionSvc.o : $(src)AntiNeutronCorrectionSvc.cxx
	$(cpp_echo) $(src)AntiNeutronCorrectionSvc.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(lib_AntiNeutronCorrectionSvc_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(use_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(lib_AntiNeutronCorrectionSvc_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(AntiNeutronCorrectionSvc_cxx_cppflags)  $(src)AntiNeutronCorrectionSvc.cxx
endif
endif

else
$(bin)AntiNeutronCorrectionSvc_dependencies.make : $(AntiNeutronCorrectionSvc_cxx_dependencies)

$(bin)AntiNeutronCorrectionSvc_dependencies.make : $(src)AntiNeutronCorrectionSvc.cxx

$(bin)$(binobj)AntiNeutronCorrectionSvc.o : $(AntiNeutronCorrectionSvc_cxx_dependencies)
	$(cpp_echo) $(src)AntiNeutronCorrectionSvc.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(lib_AntiNeutronCorrectionSvc_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(use_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(lib_AntiNeutronCorrectionSvc_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(AntiNeutronCorrectionSvc_cxx_cppflags)  $(src)AntiNeutronCorrectionSvc.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),AntiNeutronCorrectionSvcclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)AntiNeutronTrk.d

$(bin)$(binobj)AntiNeutronTrk.d :

$(bin)$(binobj)AntiNeutronTrk.o : $(cmt_final_setup_AntiNeutronCorrectionSvc)

$(bin)$(binobj)AntiNeutronTrk.o : $(src)AntiNeutronTrk.cxx
	$(cpp_echo) $(src)AntiNeutronTrk.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(lib_AntiNeutronCorrectionSvc_pp_cppflags) $(AntiNeutronTrk_pp_cppflags) $(use_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(lib_AntiNeutronCorrectionSvc_cppflags) $(AntiNeutronTrk_cppflags) $(AntiNeutronTrk_cxx_cppflags)  $(src)AntiNeutronTrk.cxx
endif
endif

else
$(bin)AntiNeutronCorrectionSvc_dependencies.make : $(AntiNeutronTrk_cxx_dependencies)

$(bin)AntiNeutronCorrectionSvc_dependencies.make : $(src)AntiNeutronTrk.cxx

$(bin)$(binobj)AntiNeutronTrk.o : $(AntiNeutronTrk_cxx_dependencies)
	$(cpp_echo) $(src)AntiNeutronTrk.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(lib_AntiNeutronCorrectionSvc_pp_cppflags) $(AntiNeutronTrk_pp_cppflags) $(use_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(lib_AntiNeutronCorrectionSvc_cppflags) $(AntiNeutronTrk_cppflags) $(AntiNeutronTrk_cxx_cppflags)  $(src)AntiNeutronTrk.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),AntiNeutronCorrectionSvcclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)AntiNeutronCorrectionSvc_entries.d

$(bin)$(binobj)AntiNeutronCorrectionSvc_entries.d :

$(bin)$(binobj)AntiNeutronCorrectionSvc_entries.o : $(cmt_final_setup_AntiNeutronCorrectionSvc)

$(bin)$(binobj)AntiNeutronCorrectionSvc_entries.o : $(src)components/AntiNeutronCorrectionSvc_entries.cxx
	$(cpp_echo) $(src)components/AntiNeutronCorrectionSvc_entries.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(lib_AntiNeutronCorrectionSvc_pp_cppflags) $(AntiNeutronCorrectionSvc_entries_pp_cppflags) $(use_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(lib_AntiNeutronCorrectionSvc_cppflags) $(AntiNeutronCorrectionSvc_entries_cppflags) $(AntiNeutronCorrectionSvc_entries_cxx_cppflags) -I../src/components $(src)components/AntiNeutronCorrectionSvc_entries.cxx
endif
endif

else
$(bin)AntiNeutronCorrectionSvc_dependencies.make : $(AntiNeutronCorrectionSvc_entries_cxx_dependencies)

$(bin)AntiNeutronCorrectionSvc_dependencies.make : $(src)components/AntiNeutronCorrectionSvc_entries.cxx

$(bin)$(binobj)AntiNeutronCorrectionSvc_entries.o : $(AntiNeutronCorrectionSvc_entries_cxx_dependencies)
	$(cpp_echo) $(src)components/AntiNeutronCorrectionSvc_entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(lib_AntiNeutronCorrectionSvc_pp_cppflags) $(AntiNeutronCorrectionSvc_entries_pp_cppflags) $(use_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(lib_AntiNeutronCorrectionSvc_cppflags) $(AntiNeutronCorrectionSvc_entries_cppflags) $(AntiNeutronCorrectionSvc_entries_cxx_cppflags) -I../src/components $(src)components/AntiNeutronCorrectionSvc_entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),AntiNeutronCorrectionSvcclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)AntiNeutronCorrectionSvc_load.d

$(bin)$(binobj)AntiNeutronCorrectionSvc_load.d :

$(bin)$(binobj)AntiNeutronCorrectionSvc_load.o : $(cmt_final_setup_AntiNeutronCorrectionSvc)

$(bin)$(binobj)AntiNeutronCorrectionSvc_load.o : $(src)components/AntiNeutronCorrectionSvc_load.cxx
	$(cpp_echo) $(src)components/AntiNeutronCorrectionSvc_load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(lib_AntiNeutronCorrectionSvc_pp_cppflags) $(AntiNeutronCorrectionSvc_load_pp_cppflags) $(use_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(lib_AntiNeutronCorrectionSvc_cppflags) $(AntiNeutronCorrectionSvc_load_cppflags) $(AntiNeutronCorrectionSvc_load_cxx_cppflags) -I../src/components $(src)components/AntiNeutronCorrectionSvc_load.cxx
endif
endif

else
$(bin)AntiNeutronCorrectionSvc_dependencies.make : $(AntiNeutronCorrectionSvc_load_cxx_dependencies)

$(bin)AntiNeutronCorrectionSvc_dependencies.make : $(src)components/AntiNeutronCorrectionSvc_load.cxx

$(bin)$(binobj)AntiNeutronCorrectionSvc_load.o : $(AntiNeutronCorrectionSvc_load_cxx_dependencies)
	$(cpp_echo) $(src)components/AntiNeutronCorrectionSvc_load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(AntiNeutronCorrectionSvc_pp_cppflags) $(lib_AntiNeutronCorrectionSvc_pp_cppflags) $(AntiNeutronCorrectionSvc_load_pp_cppflags) $(use_cppflags) $(AntiNeutronCorrectionSvc_cppflags) $(lib_AntiNeutronCorrectionSvc_cppflags) $(AntiNeutronCorrectionSvc_load_cppflags) $(AntiNeutronCorrectionSvc_load_cxx_cppflags) -I../src/components $(src)components/AntiNeutronCorrectionSvc_load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: AntiNeutronCorrectionSvcclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(AntiNeutronCorrectionSvc.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

AntiNeutronCorrectionSvcclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library AntiNeutronCorrectionSvc
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)AntiNeutronCorrectionSvc$(library_suffix).a $(library_prefix)AntiNeutronCorrectionSvc$(library_suffix).s? AntiNeutronCorrectionSvc.stamp AntiNeutronCorrectionSvc.shstamp
#-- end of cleanup_library ---------------
