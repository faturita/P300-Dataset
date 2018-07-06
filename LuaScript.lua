function initialize(box)

io.write("initialize has been called\n")

socket = require('socket')
print(socket._VERSION)

--client=socket.tcp()
--client:connect('www.itba.edu.ar', 80)
--cookie=client:receive()
--print(cookie)

-- change here to the host an port you want to contact
local host, port = "localhost", 7788
-- load namespace
-- convert host name to ip address
local ip = assert(socket.dns.toip(host))
-- create a new UDP object
local udp = assert(socket.udp())
-- contact daytime host
assert(udp:sendto("anything", ip, port))
-- retrieve the answer and print results

udp:settimeout(10)


io.write(assert(udp:receive()))



io.write(string.format("At time %f on input 1 got stimulation id:%s date:%s duration:%s\n", 33.2,'ddd','sss','sss'))
end


-- this function is called once by the box
function process(box)
        io.write("process has been called\n")
      count = 0

        -- enters infinite loop
        -- cpu will be released with a call to sleep
        -- at the end of the loop
        while true do

                -- gets current simulated time
                t = box:get_current_time()

                -- loops on all inputs of the box
               -- loops on every received stimulation for a given input
               for stimulation = 1, box:get_stimulation_count(1) do

                        -- logs the received stimulation
                        --io.write(string.format("At time %f on input 1 got stimulation id:%s date:%s duration:%s\n", t, get_stimulation(1, 1)))

                  --check witch stimulation it is
                  emit = 0
                  id,date,duration = box:get_stimulation(1, 1)
                  --io.write(string.format("Stimulation is : id:%d date:%f duration:%f\n", id,date,duration))

                  if (id==32769) then
                     --io.write(string.format("Stimulation is OVTK_StimulationId_ExperimentStart\n"))
                  end
                  if (id==32770) then
                     --io.write(string.format("Stimulation is OVTK_StimulationId_ExperimentStop\n"))
                  end
                  if (id==32771) then
                     --io.write(string.format("Stimulation is OVTK_StimulationId_SegmentStartOVTK_StimulationId_TrialStop\n"))
                  end
                  if (id==32772) then
                     --io.write(string.format("Stimulation is a OVTK_StimulationId_SegmentStop\n"))
                  end
                  if (id==32773) then
                     --io.write(string.format("Stimulation is a OVTK_StimulationId_TrialStart\n"))
                  end
                  if (id==32774) then
                     --io.write(string.format("Stimulation is a OVTK_StimulationId_TrialStop\n"))
                  end
                  if (id==32774) then
                     --io.write(string.format("Stimulation is a OVTK_StimulationId_TrialStop\n"))
                  end
                  if (id==32775) then
                     --io.write(string.format("Stimulation is a OVTK_StimulationId_BaselineStart\n"))
                  end
                  if (id==32776) then
                     --io.write(string.format("Stimulation is a OVTK_StimulationId_BaselineStop\n"))
                  end
                  if (id==32777) then
                     io.write(string.format("Stimulation is a OVTK_StimulationId_RestStart\n"))
                     emit = 1
                  end
                  if (id==32778) then
                     --io.write(string.format("Stimulation is a OVTK_StimulationId_RestStop\n"))
                     count=count+1
                     io.write(string.format("next Letter, idx = %d\n",count))
                  end
                  if (id==32779) then
                     --io.write(string.format("Stimulation is a Visual Start\n"))
                     --emit = 1
                  end
                  if (id==32772) then
                     --io.write(string.format("Stimulation is a Segment Stop\n"))
                     --count=count+1
                     --io.write(string.format("next Letter, idx = %d\n",count))
                  end

                        -- discards it
                        box:remove_stimulation(1, 1)

                        -- add triggers : new OVTK_StimulationId_Label_XX stimulation
                  if (emit==1) then
                     if (count < #CommandInput) then
                        io.write(string.format("case r = %s, c= %s\n",CommandInput[count+1][1][1],CommandInput[count+1][2][2]))
                        box:send_stimulation(1, CommandContextRow[CommandInput[count+1][1][1]], t, duration)
                        box:send_stimulation(1, CommandContextCol[CommandInput[count+1][2][2]], t, duration)
                     else
                        io.write("out of CommandInput\n")
                     end
                  end

                  --resend the base stimulation
                  --box:send_stimulation(1, id, t, duration)
               end

                -- releases cpu
                box:sleep()
        end
end



-- this function is called when the box is uninitialized
function uninitialize(box)
        io.write("uninitialize had been called\n")

end
