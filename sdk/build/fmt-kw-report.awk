BEGIN   {FS=";"}

        {sub(/.*\/sdk\//, "/sdk/", $1)  # remove long path from file name
         if ($13 == "Fixed") print "\nFIXED ISSUE"
         else if ($13 == "New") print "\nNEW ISSUE"
         else print "\nEXISTING ISSUE"
         print "FILE: " $1
         print "LINE: " $2
         print "TYPE: " $7
         print "ERROR: " $11
         print "URL: " $17
        }
        

