"""Use the caller's existing docker membership even if its user manager predates that membership."""
import grp,os,shlex,sys
argv=['/usr/bin/docker',*sys.argv[1:]]
if grp.getgrnam('docker').gr_gid in [os.getgid(),*os.getgroups()]:
 os.execv(argv[0],argv)
# sg validates membership through the system group database; no socket chmod or daemon changes.
os.execv('/usr/bin/sg',['sg','docker','-c',shlex.join(argv)])
