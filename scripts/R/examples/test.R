library(openbabelR)

convertTest <- function(){
	smi = "cc1cccc1c"

	inStr = openbabelR:::istreamFromString(smi)
	outStr = openbabelR:::ostreamToString()

	conv = openbabelR:::OBConversion(inStr,outStr)

	message("setting formats")
	openbabelR:::OBConversion_SetInAndOutFormats(conv,"SMI","SDF")

	message("converting")
	openbabelR:::OBConversion_Convert(conv)

	str = openbabelR:::stringFromOstream(outStr)
	print(str)

}
descriptorTest <- function(){
	mol = openbabelR:::OBMol()

	smi = "cc1cccc1c" 
	inStr = openbabelR:::istreamFromString(smi)

	conv = openbabelR:::OBConversion(inStr)

	message("setting formats")
	openbabelR:::OBConversion_SetInAndOutFormats(conv,"SMI","SDF")

	message("reading mol")
	openbabelR:::OBConversion_Read(conv,mol)
	message("finding type")
	desc = openbabelR:::OBDescriptor_FindType("logP")

	message("predicting")
	p= openbabelR:::OBDescriptor_Predict(desc,mol)
	message("value: ",p)

	strDesc = openbabelR:::OBDescriptor_FindType("cansmi")
	cansmi = openbabelR:::stringp()
	openbabelR:::OBDescriptor_GetStringValue(strDesc,mol,cansmi$cast())
	message("cansmi: ",cansmi$value())

}

convertTest()
descriptorTest()
