py.install_sources(
    '__init__.py',
    'density.py',
    'weighting.py',
    'zonalize.py',
    subdir: 'rwguide/zonalization'
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
    subdir: 'rwguide/zonalization',
    install: true,
    dependencies: [_ext_dep, py.dependency()]
)

