# encoding: utf-8
'''
A tool for automatic transfer of minimal relevant data from the batch server to
local storage.

Created on 15 Jan 2019

@author: Eugen Rožić
'''
import os
import shutil
import time
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('user_host', type=str,
                    help='user@host for SCP')
parser.add_argument('-s', '--search', type=str, default=r'growth',
                    help='the root search dir')
args = parser.parse_args()

scp_prefix = args.user_host
user, host = scp_prefix.split('@')
remote_local = r'/run/user/1000/gvfs/sftp:host={:s},user={:s}/'.format(host, user)
if host.startswith('myriad'):
    remote_homedir = remote_local + 'lustre/home/ucapero'
    remote_rel_basedir = r'Scratch/amyloids'
elif host.startwith('lemon'):
    remote_homedir = remote_local + 'storage/users/erozic'
    remote_rel_basedir = r'simulations/amyloids/data'
else:
    raise Exception('Unsupported host ({:s})!'.format(host))

remote_basedir = os.path.join(remote_homedir, remote_rel_basedir)
local_basedir = r'/media/data_ntfs/PhD/simulation stuff/amyloids/data'

def remote_to_local(path):
    return os.path.join(local_basedir, os.path.relpath(path, remote_basedir))

def remote_to_scp(path):
    return r"{:s}:'{:s}'".format(scp_prefix, os.path.relpath(path, remote_homedir))

def qesc(str):
    if not '"' in str:
        return r'"{:s}"'.format(str)
    elif not "'" in str:
        return r"'{:s}'".format(str)
    else:
        raise Exception("Can't escape a string that has both single and double quotation marks!")

rootsearchdir = args.search
rootnode = os.path.join(remote_basedir, rootsearchdir)

start = time.time()
for node, dirs, files in os.walk(rootnode):
    if os.path.basename(node) == 'old':
        del dirs
        continue
    for filename in files:
        basename, ext = os.path.splitext(filename)
        if not ext == '.cfg':
            continue
        if os.path.basename(node) == basename:
            # a dir containing a .cfg with the same name is a cue for it being
            # processed and to be skipped (shouldn't contain nested relevant dirs)
            del dirs
            break
        if basename in dirs:
            dirs.remove(basename) #no need to traverse this folder later
            shutil.move(os.path.join(node, filename),
                        os.path.join(node, basename, filename))
            srcdir = os.path.join(node, basename)
            print "processing", os.path.relpath(srcdir, remote_basedir), "...",
            destdir = remote_to_local(srcdir)
            if not os.path.exists(destdir):
                os.makedirs(destdir)
            scp_srcdir = remote_to_scp(srcdir)
            os.system(r"scp {:s} {:s}".format(qesc(os.path.join(scp_srcdir, '*.cfg')),
                                              qesc(destdir)))
            os.system(r"scp {:s} {:s}".format(qesc(os.path.join(scp_srcdir, '*_last_dump')),
                                              qesc(destdir)))
            os.system(r"scp {:s} {:s}".format(qesc(os.path.join(scp_srcdir, '*_cluster_data')),
                                              qesc(destdir)))
            print "done!"
secs = time.time() - start
hrs = int(secs/3600); secs -= hrs*3600
mins = int(secs/60); secs -= mins*60
print "time: {:02d}:{:02d}:{:05.2f}".format(hrs, mins, secs)
