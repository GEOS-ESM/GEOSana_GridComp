#!/usr/bin/env tcsh

set echo
set debug = ""
foreach arg ($argv)
    if ("$arg" == "-db") then
       set debug = $arg
    endif
end

@ stat = 0
foreach file (test*.py)
    set echo
    $file -v $debug
    @ stat += $status
    unset echo
    sleep 2
end

if ($stat > 0) then
    set msg = "FAILURE"
else
    set msg = "SUCCESS"
endif

echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "status = $stat"
echo $msg
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo
