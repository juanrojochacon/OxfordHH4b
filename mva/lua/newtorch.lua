local nn = require("nn")
local data = require('data')
local gnuplot = require('gnuplot')

local targetfile = arg[1]

print("Reading data ".. targetfile)
local totalset = data.readfile(targetfile, 0.1)
print(collectgarbage("count")/1024, "MB used")
data.whiten(totalset)
data.aliastable(totalset)

-- Generate dataset in torch trainer format
local torchset = data.torch(totalset,1E4)

-- Specify neural network architecture
mlp = nn.Sequential()
mlp:add( nn.Linear(totalset.nKin, 10) ) -- 22 input, 10 hidden units
mlp:add( nn.Sigmoid() )
mlp:add( nn.Linear(10, 1) ) 
mlp:add( nn.Sigmoid() ) 

-- Use cross-entropy as cost function for SGD
criterion = nn.BCECriterion()
trainer = nn.StochasticGradient(mlp, criterion)
trainer.learningRate = 0.01
trainer.maxIteration = 25
trainer:train(torchset)

-- Write output to file
local outfile = "fit_new_" .. targetfile
io.output(outfile)
for i=1,totalset.nDat,1 do
	io.write(totalset.dataSource[i],' ', totalset.dataOutput[i],' ', totalset.dataWeight[i],' ', mlp:forward(totalset.dataInputs[i])[1],'\n')
end


-------------------------------------------------------------------------------------------------------
-- Compute Signal/Background plots

local lumi = 3000 -- HL-LHC luminosity

local outputs = mlp:forward(totalset.dataInputs)
local function passWeight(coords, mask)
	local out = outputs[mask]
	local wgt = totalset.dataWeight[mask]
	return coords:apply( function(v) return lumi*torch.sum(wgt[out:gt(v)]) end ) 
end

local coords = torch.linspace(0,1,50)
local sigThr = passWeight(coords:clone(), totalset.dataOutput)
local bkgThr = passWeight(coords:clone(), totalset.dataOutput:ne(1))

local SB  = torch.cdiv(sigThr, bkgThr)
local SSB = torch.cdiv(sigThr, torch.sqrt(bkgThr))

gnuplot.epsfigure('sb'..targetfile..'.eps')
gnuplot.plot({coords,SB})
gnuplot.plotflush()

gnuplot.epsfigure('ssb'..targetfile..'.eps')
gnuplot.plot({coords,SSB})
gnuplot.plotflush()

torch.save("package.dat", totalset)