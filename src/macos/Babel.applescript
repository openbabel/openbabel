-- Babel.applescript

-----------------------------------
-- Open Babel Mac OS X GUI, version 1.99.1
-- Copyright © 2002 Geoffrey R. Hutchison
-- This program (application and code) is free software;
-- you can redistribute it and/or modify it under the terms of the 
-- GNU General Public License as published by
-- the Free Software Foundation version 2 of the License.
-- See <http://www.gnu.org/licenses/licenses.html#GPL>

-- For more information on Open Babel, see <http://openbabel.sourceforge.net/>
------------------------------------

on launched theObject
	set visible of progress indicator "progress" of window "main" to false
	--	set visible of window "main" to true
	
	set inputMenu to menu "inputMenu" of box "inputBox" of window "main"
	delete every menu item of menu inputMenu
	make new menu item at beginning of menu items of menu inputMenu with properties {title:"Automatically Detect from Extension"}
	make new menu item at end of menu items of menu inputMenu with properties {title:"Hello World!"}
end launched

on clicked theObject
	if theObject is equal to button "openButton" of box "inputBox" of window "main" then
		doOpen("")
	else if theObject is equal to button "saveButton" of box "outputBox" of window "main" then
		doSave("untitled")
	else if theObject is equal to button "addH" of box "optionsBox" of window "main" then
		if (state of button "addH" of box "optionsBox" of window "main") is 1 then
			set enabled of button "adjustPH" of box "optionsBox" of window "main" to true
			set enabled of text field "pH" of box "OptionsBox" of window "main" to true
			set contents of text field "pH" of box "OptionsBox" of window "main" to "7.0"
		else
			set enabled of button "adjustPH" of box "optionsBox" of window "main" to false
			set enabled of text field "pH" of box "optionsBox" of window "main" to false
		end if
	else if theObject is equal to button "CancelButton" of window "main" then
		quit
	end if
end clicked

on choose menu item theObject
	set menuItemTitle to title of theObject as string
	if menuItemTitle = "Open..." then
		doOpen("")
	else if menuItemTitle = "Save As..." then
		doSave("untitled")
	end if
end choose menu item

on should quit after last window closed theObject
	return true
end should quit after last window closed

on doOpen(fileName)
	if the contents of fileName is equal to "" then
		set macFile to choose file with prompt "Open Babel Input file"
	else
		set macFile to fileName
	end if
	set contents of text field "inputPath" of box "inputBox" of window "main" to macFile
	if contents of text field "inputPath" of box "inputBox" of window "main" is not equal to "" then
		set enabled of text field "outputPath" of box "outputBox" of window "main" to true
		set enabled of button "saveButton" of box "outputBox" of window "main" to true
		tell application "Finder"
			set macDir to container of macFile
		end tell
		set unixFile to POSIX path of (macFile as string)
		set unixDir to POSIX path of (macDir as string)
	end if
end doOpen

on doSave(fileName)
	set saveFile to choose file name with prompt "Save Open Babel File to" default name fileName
	set contents of text field "outputPath" of box "outputBox" of window "main" to saveFile
	if contents of text field "outputPath" of box "outputBox" of window "main" is not equal to "" then
		set enabled of button "Convert" of window "main" to true
		tell application "Finder"
			set macDir to container of macFile
		end tell
		set unixFile to POSIX path of (macFile as string)
		set unixDir to POSIX path of (macDir as string)
	end if
end doSave

