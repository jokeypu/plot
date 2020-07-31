# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.plot.Debug:
/Users/jokey/Panda/plot/xcode/Debug/plot:
	/bin/rm -f /Users/jokey/Panda/plot/xcode/Debug/plot


PostBuild.plot.Release:
/Users/jokey/Panda/plot/xcode/Release/plot:
	/bin/rm -f /Users/jokey/Panda/plot/xcode/Release/plot


PostBuild.plot.MinSizeRel:
/Users/jokey/Panda/plot/xcode/MinSizeRel/plot:
	/bin/rm -f /Users/jokey/Panda/plot/xcode/MinSizeRel/plot


PostBuild.plot.RelWithDebInfo:
/Users/jokey/Panda/plot/xcode/RelWithDebInfo/plot:
	/bin/rm -f /Users/jokey/Panda/plot/xcode/RelWithDebInfo/plot




# For each target create a dummy ruleso the target does not have to exist
