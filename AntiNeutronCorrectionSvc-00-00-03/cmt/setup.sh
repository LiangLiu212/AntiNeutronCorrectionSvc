# echo "setup AntiNeutronCorrectionSvc AntiNeutronCorrectionSvc-00-00-03 in /home/lliu/boss/Workarea-7.0.8/Analysis"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtAntiNeutronCorrectionSvctempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtAntiNeutronCorrectionSvctempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=AntiNeutronCorrectionSvc -version=AntiNeutronCorrectionSvc-00-00-03 -path=/home/lliu/boss/Workarea-7.0.8/Analysis  -no_cleanup $* >${cmtAntiNeutronCorrectionSvctempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=AntiNeutronCorrectionSvc -version=AntiNeutronCorrectionSvc-00-00-03 -path=/home/lliu/boss/Workarea-7.0.8/Analysis  -no_cleanup $* >${cmtAntiNeutronCorrectionSvctempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtAntiNeutronCorrectionSvctempfile}
  unset cmtAntiNeutronCorrectionSvctempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtAntiNeutronCorrectionSvctempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtAntiNeutronCorrectionSvctempfile}
unset cmtAntiNeutronCorrectionSvctempfile
return $cmtsetupstatus

