$originDir = Get-Location;
$srcDir = $originDir.path + "\VDRCQtLib\src"
$win32Dir = $originDir.path + "\Win32"
$x64Dir = $originDir.path + "\x64"

[Environment]::SetEnvironmentVariable('VDRCQTDIR', $originDir, 'User')
[Environment]::SetEnvironmentVariable('VDRCQTSRCDIR', $srcDir, 'User')
[Environment]::SetEnvironmentVariable('VDRCQTWIN32DIR', $win32Dir, 'User')
[Environment]::SetEnvironmentVariable('VDRCQTx64DIR', $x64Dir, 'User')
