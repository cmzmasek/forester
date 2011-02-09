package queue;

# Process queueing
# SRE, Wed Sep  2 14:37:14 1998
# CVS $Id: queue.pm,v 1.1.1.1 2005/03/22 08:35:51 cmzmasek Exp $
# Master copy: see src/queue (CVS controlled, separate from pfamserver)
#
# Written for Pfam web server; suited for queuing any set of commands.
#
# API:
# 
# $mypid = $$;
# $delay_in_seconds = 2;
#
# $nqueued = &queue::CheckQueue("pfamqueue", "username", "/tmp");
# print ("There are $nqueued jobs ahead of you in line\n");
# &queue::WaitInQueue("pfamqueue", "username", "/tmp", $mypid, $delay_in_seconds);
# print ("Our turn! Working...\n");
# (do stuff)
# &queue::RemoveFromQueue("pfamqueue", "username", "/tmp", $mypid);
#
# queuedir is a directory where the script has write permissions;
# typically a tmp directory of some sort.
#


################################################################
# PFAMSERVER - The Washington University/St. Louis Pfam web server
# Copyright (C) 1995-1999 Washington University School of Medicine
# Copyright (C) 1995-1999 Sanger Centre/Genome Research Ltd.
# Copyright (C) 1998-1999 Karolinska Institutet Center for Genomics Research
# All Rights Reserved
# 
#     This source code is distributed under the terms of the
#     GNU General Public License. See the files COPYRIGHT and LICENSE
#     for details.
# 
################################################################
# RCS $Id: queue.pm,v 1.1.1.1 2005/03/22 08:35:51 cmzmasek Exp $


# WaitInQueue() - add a process id to a queue, wait for turn
#
# Arguments: queue    - name of queue (prefix of queue stamp)
#            username - name of user (middle part of queue stamp)
#            queuedir - directory to keep queue stamps in
#            mypid    - our process id 
#            delay    - number of seconds between checking queue status
#        
# Note: When it checks the queue, if a stamp is present that 
#       doesn't seem to correspond to a running process (ps -a),
#       it deletes the stamp. This protects against crashed processes
#       freezing all subsequent jobs.
#
# example: &WaitInQueue("pfamqueue", "/tmp", $mypid, 2);
#
# Returns 1 on success, 0 on failure.
#
# NOTE: You may have to set the ps command in WaitInQueue.
#       It must return all running processes.
#
sub WaitInQueue
{
    local($queue, $username, $queuedir, $mypid, $delay) = @_;
    local(@newqueue, @queuelist, %mark);
    local(*STAMP, *QUEUEDIR);
    local(%is_running);
    local(@output, $line, $pid, $waiting);

				# get list of other guys who are working
    opendir(QUEUEDIR, $queuedir);
    @queuelist = grep(/$queue\.\S*\.\d+/, readdir(QUEUEDIR));
    closedir(QUEUEDIR);
				# make stamp for our pid
    if ($username eq "") { $username = "unknown"; }
    open(STAMP, ">$queuedir/$queue.$username.$mypid") || return 0;
    close(STAMP);
				# wait for our turn
    while (1) 
    {
	if ($#queuelist == -1) { last; } # nobody ahead of us; our turn!
	sleep($delay);
				# get list of running processes
	%is_running = 0;
	@output = split(/^/, `ps -ax`);
	foreach $line (@output) {
	    $line =~ /\s*(\d+)/;
	    $is_running{$1} = 1;
	}
				# verify that the guys we're waiting for
				# are still running, and haven't crashed.
				# if they have, reap their stamps, and their
	                        # tmp files.
	foreach $waiting (@queuelist) {
	    ($name, $pid) = ($waiting =~ /$queue\.(\S*)\.(\d+)/);
	    if (! $is_running{$pid}) { unlink "$queuedir/$queue.$name.$pid"; }
	}

	# get new list of queued jobs ahead of us.
        # ignore guys who came in after we  grabbed our initial queue list; 
        # they're waiting for *us*. The crazed greps are the Perl-y
	# way of computing an intersection between two arrays.
	#
	opendir(QUEUEDIR, $queuedir);
	@newqueue = grep(/$queue\.\S*\.\d+/, readdir(QUEUEDIR));
	closedir(QUEUEDIR);
	%mark = 0;
	grep($mark{$_}++,@queuelist);
	@queuelist = grep($mark{$_},@newqueue);
    }

    1;				# time to run!
}


# CheckQueue() - return total number of processes working, other than us
#                and the total that this particular username is running.
#
# Arguments: queue, username, queuedir
#
sub CheckQueue
{
    local($queue, $username, $queuedir) = @_;
    local(*QUEUEDIR, @allqueue, $nall, $nuser);

    opendir(QUEUEDIR, $queuedir);
    @allqueue = grep(/$queue\.\S*\.\d+/, readdir(QUEUEDIR));
    closedir(QUEUEDIR);

    if ($username eq "") {$username = "unknown"; }
    $nall = $nuser = 0;
    foreach $waiting (@allqueue) {
	($name, $pid) = ($waiting =~ /$queue\.(\S*)\.(\d+)/);
	$nall++;
	if ($name eq $username) { $nuser++; }
    }
    return ($nall, $nuser);
}
    

# RemoveFromQueue() - remove a pid from a queue
#
sub RemoveFromQueue
{
    local($queue, $username, $queuedir, $pid) = @_;
    if ($username eq "") {$username = "unknown"; }
    unlink "$queuedir/$queue.$username.$pid";
}

1;
