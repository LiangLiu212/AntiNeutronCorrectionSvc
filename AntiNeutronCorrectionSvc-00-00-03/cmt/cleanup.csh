# echo "cleanup AntiNeutronCorrectionSvc AntiNeutronCorrectionSvc-00-00-03 in /home/lliu/boss/Workarea-7.0.8/Analysis"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtAntiNeutronCorrectionSvctempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtAntiNeutronCorrectionSvctempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=AntiNeutronCorrectionSvc -version=AntiNeutronCorrectionSvc-00-00-03 -path=/home/lliu/boss/Workarea-7.0.8/Analysis  $* >${cmtAntiNeutronCorrectionSvctempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=AntiNeutronCorrectionSvc -version=AntiNeutronCorrectionSvc-00-00-03 -path=/home/lliu/boss/Workarea-7.0.8/Analysis  $* >${cmtAntiNeutronCorrectionSvctempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtAntiNeutronCorrectionSvctempfile}
  unset cmtAntiNeutronCorrectionSvctempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtAntiNeutronCorrectionSvctempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtAntiNeutronCorrectionSvctempfile}
unset cmtAntiNeutronCorrectionSvctempfile
exit $cmtcleanupstatus

