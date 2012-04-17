#!/bin/bash
#
# Initscript for Pathway Tools software
#
# Author:	Lukas Mueller <lam87@cornell.edu> 
#               June 8, 2007
#               based on example initscript by:
#               Miquel van Smoorenburg <miquels@cistron.nl>.
#		Ian Murdock <imurdock@gnu.ai.mit.edu>.
#
# Version:	@(#)skeleton  2.85-23  28-Jul-2004  miquels@cistron.nl
#

### BEGIN INIT INFO
# Provides:          pathway_tools.sh
# Required-Start:    $all
# Required-Stop:     
# Default-Start:     2 3 4 5
# Default-Stop:      0 1 6
# Short-Description: Start Pathway Tools at boot time
# Description:       Run Pathway Tools as a webserver at bootup
### END INIT INFO


#set -e

PATH=/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin
DESC="pathway tools"
NAME=ptools
#DAEMON=/usr/sbin/$NAME
PIDFILE=/tmp/$NAME.pid
SCRIPTNAME=/etc/init.d/$NAME.sh

# Gracefully exit if the package has been removed.
#test -x $DAEMON || exit 0

# Read config file if it is present.
#if [ -r /etc/default/$NAME ]
#then
#	. /etc/default/$NAME
#fi

#
#	Function that starts the daemon/service.
#
function d_start() {
# Warn and exit if something is already listening on port 1555.                                                                            
    if [ `netstat -a | fgrep 1555 | fgrep -c LISTEN` == 1 ] 
	then
	echo "Something is already listening on port 1555.";
	echo "Please stop whatever it is, and then try $0 again."
	exit 1;
    fi

    if [ -a $PIDFILE ]
	then
	echo "The process file $PIDFILE already exists. ";
	echo "Pathway tools already seems to have been started using this method.";
	exit 1;

    fi

## Start Xvfb as X-server display #1. If it's already been started, no harm done.
# Start Xvfb as X-server display #1. If it's already been started, kill it.
# This works in conjunction with the pathway-tools script which needs to have:                                                            
    
    export DISPLAY=localhost:1;
	pkill Xvfb;
	sleep 2;
	echo "STARTING Xvfb..."
	/usr/X11R6/bin/Xvfb :1 &
# Start GNU screen (hit Ctrl-A Ctrl-D to detach). Within it, start Pathway Tools.                                                         
	echo "STARTING screen...";

	/usr/bin/screen -m -L -d /usr/local/pathway-tools/aic-export/pathway-tools/pttols/15.5/pathway-tools -www -www-publish public -org lyco -gene-link-db SGN  -no-blast -allow-webcrawlers;

	echo "Determine PID after wait...";
	sleep 2;

	ps -o pid -C pathway-tools-runtime --no-headers > $PIDFILE;
	
	echo `cat $PIDFILE`;
	#start-stop-daemon --start --quiet --pidfile $PIDFILE \
	#	--exec /root/sgn-tools/solcyc/start_pathway_tools.sh
}

#
#	Function that stops the daemon/service.
#
function d_stop() {
	#start-stop-daemon --stop --quiet --pidfile $PIDFILE \
	#	--name $NAME



	if [ -a $PIDFILE ]; 
	    then
	    PID=`cat $PIDFILE`;	
	    echo "";
	    echo "Process ID $PID. Killing...";
	    kill -9 $PID;
	    echo "Removing PID file $PID...";
	    rm $PIDFILE;
	    echo "Done."   

	
	else 
	   echo "";
	   echo "Process file $PIDFILE not found";
	   exit 1;
	fi
	

	
	

}

#
#	Function that sends a SIGHUP to the daemon/service.
#
function d_reload() {
    #start-stop-daemon --stop --quiet --pidfile $PIDFILE \
    #--name $NAME --signal 1
    d_stop;
    sleep 2;
    d_start;
}

case "$1" in
  start)
	echo -n "Starting $DESC: $NAME"
	d_start
	echo "."
	;;
  stop)
	echo -n "Stopping $DESC: $NAME"
	d_stop
	echo "."
	;;
  #reload)
	#
	#	If the daemon can reload its configuration without
	#	restarting (for example, when it is sent a SIGHUP),
	#	then implement that here.
	#
	#	If the daemon responds to changes in its config file
	#	directly anyway, make this an "exit 0".
	#
	# echo -n "Reloading $DESC configuration..."
	# d_reload
	# echo "done."
  #;;
  restart|force-reload)
	#
	#	If the "reload" option is implemented, move the "force-reload"
	#	option to the "reload" entry above. If not, "force-reload" is
	#	just the same as "restart".
	#
	echo -n "Restarting $DESC: $NAME"
	/etc/init.d/apache2 stop
	d_stop
	sleep 1
	d_start
	sleep 120
	/etc/init.d/apache2 start
	echo "."
	;;
  *)
	# echo "Usage: $SCRIPTNAME {start|stop|restart|reload|force-reload}" >&2
	echo "Usage: $SCRIPTNAME {start|stop|restart|force-reload}" >&2
	exit 1
	;;
esac

exit 0
