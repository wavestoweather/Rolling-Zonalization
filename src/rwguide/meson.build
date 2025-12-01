py.install_sources(
    '__init__.py',
    'bindings.py',
    'xarray/__init__.py',
    'xarray/common.py',
    'xarray/hovmoeller.py',
    'xarray/pvgradient.py',
    'xarray/waveactivity.py',
    'xarray/wavenumber.py',
    'xarray/zonalization.py',
    subdir: 'rwguide',
    preserve_path: true
)

pyext = run_command('python3-config', '--extension-suffix', check: true).stdout().strip()

subdir('hovmoeller')
subdir('pvgradient')
subdir('waveactivity')
subdir('wavenumber')
subdir('zonalization')
