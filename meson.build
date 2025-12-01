project('rwguide')

add_languages('c')

py = import('python').find_installation(pure: false)

pyext = run_command('python3-config', '--extension-suffix', check: true).stdout().strip()

subdir('src/rwguide')

