#!/bin/csh

# Author James Dempsey
# Date 27 Jun 2016

# Flag the invalid portions fo the data

#echo "##--## Flag across all data ##--##"
#uvflag vis=magmo-2010-12-11.uv flagval=f options=brief "select=ant(4)"
#uvflag vis=magmo-2010-12-11.uv flagval=f options=brief "select=ant(6),pol(xx)"
echo "##--## Split by source ##--##"
uvsplit vis=magmo-2011-09-10.uv
echo "##--## Flag for specific sources ##--##"
uvflag vis=1934-638.1420 flagval=f options=brief line=channel,10,1020,1,1
uvflag vis=1934-638.1420 flagval=f options=brief line=channel,110,2430,1,1


#uvflag vis=281.710-1.104.1420 flagval=f options=brief "select=ant(4)(5)"

#uvflag vis=284.352-0.419.1420 flagval=f options=brief "select=ant(4)(3)"
#uvflag vis=284.352-0.419.1420 flagval=f options=brief "select=ant(4)(5)"

#uvflag vis=284.694-0.361.1420 flagval=f options=brief "select=ant(4)(5)"

#uvflag vis=285.337-0.002.1420 flagval=f options=brief "select=ant(4)(5)"

#uvflag vis=287.371+0.644.1420 flagval=f options=brief "select=ant(4)(5),time(01:00:00,02:00:00)"

#uvflag vis=291.270-0.719.1420 flagval=f options=brief "select=ant(4)(5)"
#uvflag vis=291.270-0.719.1420 flagval=f options=brief "select=ant(3)(4)"
#uvflag vis=291.270-0.719.1420 flagval=f options=brief "select=ant(3)(5)"
#uvflag vis=291.270-0.719.1420 flagval=f options=brief "select=ant(1)(2)"
