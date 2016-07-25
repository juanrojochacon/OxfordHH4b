--- Neural network data sources
local data = {}
local thispath = select('1', ...):match(".+%.") or ""

local utils = require('pl.utils')
require('unsup')

local function countLines(targetfile)
	local iline = 0
	for line in io.lines(targetfile) do iline = iline+1 end
	return iline
end

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

function data.whiten(set)
	print("Whitening data")
	set.dataInputs = unsup.zca_whiten(set.dataInputs)
end

-- Walker alias
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

return data