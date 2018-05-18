# -*- mode: python -*-
import sys
import os
sys.setrecursionlimit(5000)

block_cipher = None

options = [ ('v', None, 'OPTION') ]

a = Analysis(['GUI.py'],
             pathex=['C:\\users\\ncautaer\\desktop\\Python Experience\\TEM calcs v2'],
             binaries=[],
             datas = [],
             #datas=[('\Images\*.png', 'Images')],
             hiddenimports=['scipy._lib.messagestream'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='GUI',
          debug=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='GUI')
