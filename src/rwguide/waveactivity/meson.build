py.install_sources(
    '__init__.py',
    'lwa.py',
    subdir: 'rwguide/waveactivity'
)


_ext_src = custom_target(
    command: ['python3', '@INPUT0@', '--emit-c-code', '@OUTPUT@'],
    input: ['build.py', 'ext.c'],
    output: '_ext.c',
    install: false
)

_ext_dep = declare_dependency(sources: ['ext.c'])

py.extension_module(
    '_ext',
    _ext_src,
    subdir: 'rwguide/waveactivity',
    install: true,
    dependencies: [_ext_dep, py.dependency()]
)

