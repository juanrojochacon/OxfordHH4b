--- Neural network data sources
local data = {}
local thispath = select('1', ...):match(".+%.") or ""

local nn = require("nn")
local gnuplot = require('gnuplot')
require('unsup')

function data.set(nkin)
	local newset = {}
	newset._nkin = nkin
	newset._points = {}
	return newset
end

function data.readfile(targetfile)
	local function split(string)
		local tokens = {}
		for i in string.gmatch(string, "%S+") do table.insert(tokens,i) end
		return tokens
	end
	local iline = 1 
	local set = nil
	for line in io.lines(targetfile) do 
	    local tokens = split(line)
	    if iline == 1 then
	    	set = data.set(#tokens - 4)
	    else
		    local kinematics = {} for ikin=4, #line,1 do table.insert(kinematics, tonumber(tokens[ikin])) end
		    data.addpoint(set, tonumber(tokens[1]), tokens[2], tonumber(tokens[3]), kinematics)
	 	end
	 	iline = iline+1
	 end
	 return set
end

function data.addpoint(set, sig, src, wgt, kin)
	assert(#kin == set._nkin, "kinematics mismatch")
	table.insert(set._points, {signal=sig, source=src, weight=wgt, kinematics=kin} )
end

function data.normalise(set)
	print("Normalising data")
	for ikin=1,set._nkin,1 do
		local kinmax = 0
		local kinmin = math.huge
		for _,point in ipairs(set._points) do
			kinmax = math.max(kinmax,point.kinematics[ikin])
			kinmin = math.min(kinmin,point.kinematics[ikin])
		end
		for _,point in ipairs(set._points) do
			norm = (2.0*math.sqrt(3.0)/(kinmax - kinmin))
			shift = - ( norm*kinmin + math.sqrt(3) )
			point.kinematics[ikin] = point.kinematics[ikin]*norm + shift
		end	
	end
end

function data.whiten(set)
	print("Whitening data")
	local dataTensor = torch.Tensor(#set._points, set._nkin):zero()
	for k,p in ipairs(set._points) do
		for i=1,set._nkin do
			dataTensor[k][i] = p.kinematics[i]
		end
	end

	local auxTensor, _, _ = unsup.zca_whiten(dataTensor)
	print("Transform calculated")
	for k,p in ipairs(set._points) do
		for i=1,set._nkin do
			p.kinematics[i] = auxTensor[k][i]
			-- for j=1,set._nkin do
			-- 	p.kinematics[i] = p.kinematics[i] + tr[i][j]*dataTensor[k][j]
			-- end
		end
	end


	-- pT1 = dataTensor:select(2, 1)
	-- gnuplot.hist(pT1)

	-- pT2 = auxdata:select(2, 1)
	-- gnuplot.hist(pT2)


	-- gnuplot.imagesc(covMat)
   --gnuplot.plotflush()
end

function data.print(set)
	for i = 1, #set._points,1 do
		print(set._points[i].source, set._points[i].signal, set._points[i].weight)
	end
end

function data.shuffle(set)
	for i = 1, #set,1 do
		local j = math.random(#set)
		set[i], set[j] = set[j], set[i]
	end
end

function data.torch(set, lumi, overweight)
	torchset={};
	function torchset:size() return #self end
	for k,pt in ipairs(set._points) do 
		print(k)
	  local input = torch.Tensor(pt.kinematics);     -- normally distributed example in 2d
	  local output = torch.Tensor({pt.signal});
	  local eventCount = math.floor(pt.weight*lumi*(overweight^pt.signal) + 0.5)
	  -- if pt.signal == 1 then print(eventCount, pt.weight) end
	  for i=1,eventCount,1 do
		  table.insert(torchset,{input, output})
	  end
	end
	return torchset
end

return data