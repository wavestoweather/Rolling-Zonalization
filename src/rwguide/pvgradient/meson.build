py.install_sources(
    '__init__.py',
    'processing.py',
    'variables.py',
    subdir: 'rwguide/pvgradient'
)


_ext_src = custom_target(
    command: ['python3', '@INPUT0@', '--emit-c-code', '@OUTPUT0@'],
    input: ['build.py', 'ext.c'],
    output: '_ext.c',
    install: false
)

_ext_dep = declare_dependency(sources: ['ext.c'])

py.extension_module(
    '_ext',
    _ext_src,
    subdir: 'rwguide/pvgradient',
    install: true,
    dependencies: [_ext_dep, py.dependency()]
)

