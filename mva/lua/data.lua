--- HH4b Data sources
-- This file contains scripts for the management of data nTuples from
-- particle physics analysis code.

local data = {}
local thispath = select('1', ...):match(".+%.") or ""

local utils = require('pl.utils')
require('unsup')

-- A (not particuarly efficient) function to count the number of lines in a file
local function countLines(targetfile)
	local iline = 0
	for line in io.lines(targetfile) do iline = iline+1 end
	return iline
end

--- readfile
-- This function takes a target filename corresponding to a standard NTuple format,
-- and returns a 'dataset' table.

-- The dataset table has several components:
-- 	set.dataInputs: a (nData,nKinematics)-sized torch.Tensor holding the input data points
-- 	set.dataWeight: a (nData)-sized torch.Tensor holding the weight for each input data point
-- 	set.dataOutput: a (nData)-sized torch ByteTensor holding the 'output' of the datapoints (i.e 1 for signal, 0 for background)
-- 	set.dataSource: a (nData)-sized list describing the source of each data point.
-- 	set.nSig:		an integer describing how many signal points are in the sample
-- 	set.nBkg:		an integer describing now many background points are in the sample
--  set.signorm:	the overall normalisation for signal events in the training cost
--  set.bkgnorm:    the overall normalisation for background events in the training cost
function data.readfile(targetfile, sigfr)
	local set = {}

	local nData = countLines(targetfile) - 1
	print(nData .. " datapoints in target file")

	local file = io.open(targetfile, "r")
	io.input(file)

	io.open(targetfile,'r')
	local line = utils.split(io.read())
	set.nKin = #line - 4
	print(set.nKin .. " kinematics in target file")

	set.dataInputs = torch.Tensor(nData, set.nKin):zero()
	set.dataWeight = torch.Tensor(nData):zero()
	set.dataOutput = torch.ByteTensor(nData):zero()
	set.dataSource = {}
	set.nSig = 0
	set.nBkg = 0

	local sigwgt = 0
	local bkgwgt = 0

	-- Read datapoints from file, assuming the layout
	-- <signal flag [0/1]> <source description string> <datapoint weight> <kin1> ... <kinN>
	for i=1,nData do
		local line = utils.split(io.read())
		for j=1, set.nKin do set.dataInputs[i][j] = tonumber(line[j+3]) end
		set.dataOutput[i] = tonumber(line[1])
		set.dataSource[i] = line[2]
		set.dataWeight[i] = tonumber(line[3])
		if set.dataOutput[i] == 1 then 
			set.nSig = set.nSig + 1 
			sigwgt = sigwgt + set.dataWeight[i]
		else
			set.nBkg = set.nBkg + 1 
			bkgwgt = bkgwgt + set.dataWeight[i]
		end
	end

	set.nDat = set.nSig + set.nBkg

	print(sigwgt .. " signal weight")
	print(bkgwgt .. " background weight")

	set.signorm = sigfr / sigwgt
	set.bkgnorm = (1.0-sigfr) / bkgwgt

	io.close(file)
	return set
end

-- Decorrelate the kinematical variables in the dataset by ZCA-whitening
-- This is a helper function which you can use to try and simplify NN training
function data.whiten(set)
	print("Whitening data")
	set.dataInputs = unsup.zca_whiten(set.dataInputs)
end

-- This function computes the alias table for Walker's alias method.
-- The alias table is then used for selecting points at random from the
-- dataset in an efficient manner (see data.wgtpoint)
function data.aliastable(set)
	print("Computing alias table")
	local overfull = {}
	local underfull = {}
	local exactfull = {}

	local U = {}
	local K = {}
	local n = set.nDat
	for i=1,n,1 do 
		if set.dataOutput[i] == 1 then
			U[i] = n*set.dataWeight[i]*set.signorm
		else
			U[i] = n*set.dataWeight[i]*set.bkgnorm
		end
		if U[i] > 1 then table.insert(overfull, i)
		elseif U[i] < 1 then table.insert(underfull, i)
		else table.insert(exactfull,i) end
	end

	while #overfull > 0 and #underfull > 0 do
		local i = overfull[#overfull]
		local j = underfull[#underfull]
		K[j] = i
		U[i] = U[i] + U[j] - 1.0
		table.remove(underfull)
		table.insert(exactfull,j)

		if U[i] == 1 then 
			table.remove(overfull)
			table.insert(exactfull, i)
		elseif U[i] < 1 then 
			table.remove(overfull)
			table.insert(underfull, i)
		end
	end

	-- Handle remaining bin
	assert(#overfull + #underfull == 1)
	if #underfull > 0 then	U[underfull[1]] = 1 end
	if #overfull > 0 then	U[overfull[1]] = 1 end

	set.U = U
	set.K = K
end

-- Function that returns a random element index of the dataset
-- according to the alias method
function data.wgtpoint(set)
	local n = #set.U
	local x = math.random()
	local i = math.floor(n*x) + 1.0
	local y = n*x + 1.0 - i
	if (y < set.U[i]) then 
		return i 
	else 
		return set.K[i] 
	end
end

-- Generate dataset for torch standard format
function data.torch(set, size)
	torchset={};
	function torchset:size() return #self end
	for i=1,size,1 do
		local point = data.wgtpoint(set)
		local inp = set.dataInputs[point]
		local oup = set.dataOutput[point]
		table.insert(torchset, {inp, torch.Tensor({oup}) } )
	end
	return torchset
end

-- Generate a sub-batch of training elements
function data.batch(set, size)
	local inputs = torch.Tensor(size, set.nKin)
	local labels = torch.DoubleTensor(size)
	for i=1,size do
		local point = data.wgtpoint(set)
		labels[i] = set.dataOutput[point]
		inputs[i]:copy(set.dataInputs[point])
	end
	return labels, inputs
end

return data