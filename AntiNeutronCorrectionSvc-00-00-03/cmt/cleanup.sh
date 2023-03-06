# echo "cleanup AntiNeutronCorrectionSvc AntiNeutronCorrectionSvc-00-00-03 in /home/lliu/boss/Workarea-7.0.8/Analysis"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtAntiNeutronCorrectionSvctempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtAntiNeutronCorrectionSvctempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=AntiNeutronCorrectionSvc -version=AntiNeutronCorrectionSvc-00-00-03 -path=/home/lliu/boss/Workarea-7.0.8/Analysis  $* >${cmtAntiNeutronCorrectionSvctempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=AntiNeutronCorrectionSvc -version=AntiNeutronCorrectionSvc-00-00-03 -path=/home/lliu/boss/Workarea-7.0.8/Analysis  $* >${cmtAntiNeutronCorrectionSvctempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtAntiNeutronCorrectionSvctempfile}
  unset cmtAntiNeutronCorrectionSvctempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtAntiNeutronCorrectionSvctempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtAntiNeutronCorrectionSvctempfile}
unset cmtAntiNeutronCorrectionSvctempfile
return $cmtcleanupstatus

