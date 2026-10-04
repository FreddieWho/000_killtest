"""Filesystem/process event dependency, with no periodic status polling."""
import ctypes,json,os,pathlib,select
def wait_receipt(path,controller_pid):
 path=pathlib.Path(path);libc=ctypes.CDLL(None,use_errno=True)
 fd_pid=libc.syscall(434,controller_pid,0)
 def read():
  try:return json.loads(path.read_text())
  except (FileNotFoundError,json.JSONDecodeError):return None
 try:
  while True:
   result=read()
   if result is not None:
    if result.get('exit_code')!=0:raise RuntimeError('Dependency failed: '+str(path))
    return result
   if fd_pid<0:raise RuntimeError('Controller absent and dependency incomplete: '+str(path))
   parent=path.parent
   while not parent.exists():parent=parent.parent
   fd=libc.inotify_init1(os.O_CLOEXEC)
   if fd<0:raise OSError(ctypes.get_errno(),'inotify')
   try:
    if libc.inotify_add_watch(fd,os.fsencode(parent),0x8|0x80|0x100)<0:raise OSError(ctypes.get_errno(),'watch')
    if read() is not None:continue
    ready=select.select([fd,fd_pid],[],[])[0]
    if fd_pid in ready and read() is None:raise RuntimeError('Controller exited before dependency completed: '+str(path))
    if fd in ready:os.read(fd,65536)
   finally:os.close(fd)
 finally:
  if fd_pid>=0:os.close(fd_pid)
