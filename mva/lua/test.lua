local data = require('data')
local ann = require('ann')

local totalset = data.set(22)
local targetfile = "bstNTuple.dat"

print("Reading data")
data.readfile(totalset, targetfile)
print(#totalset._points .. " points read")

print("Normalising data")
data.normalise(totalset)
data.shuffle(totalset)

print("Initialising network")
local network =  ann.new( {22,5,1} )

for i=1,100,1 do
	local totCE = 0
	for _,point in ipairs(totalset._points) do
		local output = ann.compute(network, point.kinematics)
		local t   = point.signal
		local tpr = unpack(output)
		local wgt = point.weight
		local ce  = -(t*math.log(tpr)+(1.0-t)*math.log(1.0-tpr))
		local dE  = {-t/tpr + (1.0-t)/(1.0-tpr)}
		totCE = totCE + ce

		ann.backprop(network, dE, 0.01)
	end
	print(i,totCE)
end

