import commands
import sys
import glob

wd = sys.argv[1]
things = glob.glob(wd + '/*')

for thing in things:
    if len(thing.split('_')) == 7:
        thing_split = thing.split('_')
        new_thing = '_'.join(thing_split[:6]) + '_wedge_output'
        print 'changing ' + thing + ' to ' + new_thing
        commands.getoutput('mv %s %s' % (thing, new_thing))

    else:
        pass

