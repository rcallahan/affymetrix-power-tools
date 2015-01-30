#
# affy/sdk/build/cppunit_funcs.sh ---
#
# $Id$
#

# cppunit_rm_file LIST_OF_FILENAMES
function cppunit_rm_file () {
  # avoid the error reports from rm.
  for file in "${@}"
  do
    if [ -e "${file}" ]
    then      
      rm "${file}"
    fi
  done
}

# cppunit_set_outputfile CPPUNIT_XMLFILENAME
# sets the final output filename and inits counters.
function cppunit_set_outputfile () {
  # remember for later and gen tmp filenames
  cppunit_output_file="$1"
  cppunit_tmp_file="$1.tmp-$$"
  cppunit_tmp_failed="${cppunit_tmp_file}.failed"
  cppunit_tmp_passed="${cppunit_tmp_file}.passed"
  # 
  cppunit_total_tests=0
  cppunit_total_failed=0
  cppunit_total_passed=0
  # this has the name of the active test.
  cppunit_test_name=""
  #
  cppunit_rm_file "${cppunit_output_file}" "${cppunit_tmp_failed}"  "${cppunit_tmp_passed}"
  # cause we will cat these later.
  touch "${cppunit_tmp_failed}" "${cppunit_tmp_passed}"
  #
  cppunit_xsl_file=${sdk_root}/build/cppunit2junit.xsl
}

# cppunit_test_start TESTNAME FILE? LINE?
# FILE and LINE are optional
# Sets the name of this test and where it was defined.
function cppunit_test_start () {
  if [ "${cppunit_test_name}" != "" ]
   then 
    echo "cppunit_test_start(): No result for the last test '${cppunit_test_name}'!..."
  fi  
  #
  if [ "${1}" == "" ]
  then 
    echo "cppunit_test_start(): error: no name for test."
    cppunit_test_name="${BASH_SOURCE[1]}-${cppunit_total_tests}"
  else
    cppunit_test_name="$1"
  fi
  # @todo Check for "::" in the name (For XSL conversion with ${cppunit_xsl_file})
  #
  # the top of the stack is 0.
  # echo "BASH_SOURCE=${BASH_SOURCE[*]}"
  # echo "BASH_LINENO=${BASH_LINENO[*]}"
  #
  if [ "${2}" != "" ]
  then
    cppunit_test_file="$2"
  else
    cppunit_test_file="${BASH_SOURCE[0]}"
  fi
  if [ "$3" != "" ]
  then
    cppunit_test_line="$3"
  else
    cppunit_test_line="${BASH_LINENO[0]}"
  fi
}

# append the passed test to the file and bump the counters
# no args are needed, cause its ok.
function cppunit_test_passed () {
  #
  cppunit_total_tests=$(( ${cppunit_total_tests} + 1 ))
  cppunit_total_passed=$(( ${cppunit_total_passed} + 1 ))
  #
  cat >>"${cppunit_tmp_passed}" <<__EOF
  <Test id="${cppunit_total_tests}">
    <Name>${cppunit_test_name}</Name>
  </Test>
__EOF
  # let cppunit_test_start know a result was posted.
  cppunit_test_name=""
}

# cppunit_test_failed CODE MESSAGE
# CODE = a CPPunit failure code
# MESSAGE = text
# Marks this test as having failed.
function cppunit_test_failed () {
  #
  cppunit_total_tests=$(( ${cppunit_total_tests} + 1 ))
  cppunit_total_failed=$(( ${cppunit_total_failed} + 1 ))
  #
  if [ "${1}" == "" ]
  then        
    echo "cppunit_test_failed(): No failure type set"
    cppunit_test_failtype="UnsetFailure"
  else
    cppunit_test_failtype="${1}"
  fi
    
  #
  cat >>"${cppunit_tmp_failed}" <<__EOF
  <FailedTest id="${cppunit_total_tests}">
    <Name>${cppunit_test_name}</Name>
    <FailureType>${cppunit_test_failtype}</FailureType>
    <Location>
      <File>${cppunit_test_file}</File>
      <Line>${cppunit_test_line}</Line>
    </Location>
    <Message>${2}</Message>
  </FailedTest>
__EOF
  # let cppunit_test_start know a result was posted.
  cppunit_test_name=""
}

# assemble the output and clean up our tmp files.
function cppunit_finish () {
  # we now use the base tmp file name to assemble the final xml
  cppunit_rm_file "${cppunit_tmp_file}"
  #
  cat >>"${cppunit_tmp_file}" <<__EOF
<?xml version="1.0" encoding='ISO-8859-1' standalone='yes' ?>
<TestRun>
<FailedTests>
__EOF
  #
  cat "${cppunit_tmp_failed}" >>"${cppunit_tmp_file}"
  #
  cat >>"${cppunit_tmp_file}" <<__EOF
</FailedTests>
<SuccessfulTests>
__EOF
  #
  cat "${cppunit_tmp_passed}" >>"${cppunit_tmp_file}"
  #
  cat >>"${cppunit_tmp_file}" <<__EOF
</SuccessfulTests>
<Statistics>
  <Tests>${cppunit_total_tests}</Tests>
  <FailuresTotal>${cppunit_total_failed}</FailuresTotal>
  <Failures>${cppunit_total_failed}</Failures>
  <Errors>${cppunit_total_errors}</Errors>
</Statistics>
</TestRun>
__EOF
  # now that it is assembled, put it where it belongs and clean up.
  cat "${cppunit_tmp_file}" > "${cppunit_output_file}"
  cppunit_rm_file "${cppunit_tmp_file}" "${cppunit_tmp_failed}" "${cppunit_tmp_passed}"
  # Now convert the result like the SDK makefile does.
  if [ -e ${cppunit_xsl_file} ]
  then
    xsltproc ${cppunit_xsl_file} ${cppunit_output_file} > JUnitTestResults.xml
  fi
}
