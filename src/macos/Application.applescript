-- Application.applescript

-----
-- Open Babel Mac OS X GUI, version 1.99.1
-- Copyright © 2002 Geoffrey R. Hutchison
-- This program (application and code) is free software;
-- you can redistribute it and/or modify it under the terms of the 
-- GNU General Public License as published by
-- the Free Software Foundation version 2 of the License.
-- See <http://www.gnu.org/licenses/licenses.html#GPL>

-- For more information on Open Babel, see <http://openbabel.sourceforge.net/>
------

on launched theObject
	-- to be added shortly
end launched

on clicked theObject
	if theObject is equal to button "OpenButton" of box "InputBox" of window "Open Babel" then
		set macFile to choose file with prompt "Open Babel Input file"
		set contents of text field "InputPath" of box "InputBox" of window "Open Babel" to macFile
		if contents of text field "InputPath" of box "InputBox" of window "Open Babel" is not equal to "" then
			set enabled of button "Convert" of window "Open Babel" to true
			set enabled of text field "OutputPath" of box "OutputBox" of window "Open Babel" to true
			set enabled of button "SaveButton" of box "OutputBox" of window "Open Babel" to true
		end if
		tell application "Finder"
			set macDir to container of macFile
		end tell
		set unixFile to POSIX path of (macFile as string)
		set unixDir to POSIX path of (macDir as string)
	else if theObject is equal to button "AddH" of box "OptionsBox" of window "Open Babel" then
		if (state of button "AddH" of box "OptionsBox" of window "Open Babel") is 1 then
			set enabled of button "AdjustPH" of box "OptionsBox" of window "Open Babel" to true
			set enabled of text field "pH" of box "OptionsBox" of window "Open Babel" to true
			set contents of text field "pH" of box "OptionsBox" of window "Open Babel" to "7.0"
		else
			set enabled of button "AdjustPH" of box "OptionsBox" of window "Open Babel" to false
			set enabled of text field "pH" of box "OptionsBox" of window "Open Babel" to false
		end if
	else if theObject is equal to button "CancelButton" of window "Open Babel" then
		quit
	end if
end clicked

on choose menu item theObject
	set menuItemTitle to title of theObject as string
	if menuItemTitle = "Open..." then
		set macFile to choose file with prompt "Open Babel Input file"
		set contents of text field "InputPath" of box "InputBox" of window "Open Babel" to macFile
		if contents of text field "InputPath" of box "InputBox" of window "Open Babel" is not equal to "" then
			set enabled of button "Convert" of window "Open Babel" to true
			set enabled of text field "OutputPath" of box "OutputBox" of window "Open Babel" to true
			set enabled of button "SaveButton" of box "OutputBox" of window "Open Babel" to true
		end if
		tell application "Finder"
			set macDir to container of macFile
		end tell
		set unixFile to POSIX path of (macFile as string)
		set unixDir to POSIX path of (macDir as string)
	else if menuItemTitle = "Save As..." then
		set saveFile to choose file name with prompt "Save File to" default name "Testing 1"
	end if
end choose menu item

on should quit after last window closed theObject
	return true
end should quit after last window closed
