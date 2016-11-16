local nn = require("nn")
local data = require('data')
local gnuplot = require('gnuplot')
require('optim')

local targetfile = arg[1]

print("Reading data ".. targetfile)
local totalset = data.readfile(targetfile, 0.1)
print(collectgarbage("count")/1024, "MB used")

-- Setting seeds
torch.manualSeed(0)
math.randomseed(0)

--Preprocessing
data.whiten(totalset)
data.aliastable(totalset)

mlp = nn.Sequential()
mlp:add( nn.Linear(totalset.nKin, 10) ) -- 22 input, 5 hidden units
mlp:add( nn.Sigmoid() )
mlp:add( nn.Linear(10, 1) ) 
mlp:add( nn.Sigmoid() ) 

local params, gradParams = mlp:getParameters()
local criterion = nn.BCECriterion()

local optimState = {
   learningRate = 0.01,
   learningRateDecay = 0,
   weightDecay = 0,
   --momentum = 0.05
   }

local sampleSize, batchSize = 1E4,10
local ne = 200*batchSize

for epoch = 1,ne,1 do
	local sampleLabels, sampleInputs = data.batch(totalset,sampleSize)
	for ibatch=1,sampleSize/batchSize do

		local start = batchSize*(ibatch-1) + 1
		local stop  = start + (batchSize - 1)

		-- print(start, stop, batchSize, sampleSize, #sampleLabels)

		local batchLabels = sampleLabels:sub(start,stop)
		local batchInputs = sampleInputs:sub(start,stop)

		local function feval(newparams)
			if params ~= newparams then
				params:copy(newparams)
			end gradParams:zero()
			local output = mlp:forward(batchInputs)
			local erf = criterion:forward(output, batchLabels)
			mlp:backward(batchInputs, criterion:backward(output, batchLabels))
			return erf, gradParams
		end

		local _, fs = optim.nag(feval, params, optimState)
	end

	local output = mlp:forward(sampleInputs)
	local erf = criterion:forward(output, sampleLabels)
	print(epoch, erf)
end

local outfile = "fit_opt_" .. targetfile
io.output(outfile)
for i=1,totalset.nDat,1 do
	io.write(totalset.dataSource[i],' ', totalset.dataOutput[i],' ', totalset.dataWeight[i],' ', mlp:forward(totalset.dataInputs[i])[1],'\n')
end
