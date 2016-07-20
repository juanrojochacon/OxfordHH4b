local nn = require("nn")
local data = require('data')

local targetfile = arg[1]
local ikin = tonumber(arg[2])

print("Reading data ".. targetfile)

local totalset = data.readfile(targetfile)
print(#totalset._points .. " points read")
print(collectgarbage("count")/1024, "MB used")

--data.whiten(totalset)
--local torchset = data.torch(totalset, 300, 5E2)
local torchset = data.torch(totalset, 3, 5E4)

print(#torchset .. " points for training")
print(collectgarbage("count")/1024, "MB used")

data.shuffle(torchset)

mlp = nn.Sequential()
mlp:add( nn.Linear(totalset._nkin, 10) ) -- 22 input, 5 hidden units
mlp:add( nn.Sigmoid() )
mlp:add( nn.Linear(10, 1) ) 
mlp:add( nn.Sigmoid() ) 

criterion = nn.BCECriterion()
trainer = nn.StochasticGradient(mlp, criterion)
trainer.learningRate = 0.01
trainer.maxIteration = 25
trainer:train(torchset)

local outfile = "res_" .. targetfile
io.output(outfile)
for _,v in ipairs(totalset._points) do
	io.write(v.source,' ', v.signal,' ', v.weight,' ', mlp:forward(torch.Tensor(v.kinematics))[1],'\n')
end
