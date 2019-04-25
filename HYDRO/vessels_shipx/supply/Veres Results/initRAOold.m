for i=1:6

amp{i}   = reshape(vessel.forceRAO.amp{i}(:,:,1),36,19);
phase{i} = reshape(vessel.forceRAO.phase{i}(:,:,1),36,19);

end
for i=1:3

amp2{i}   = reshape(vessel.driftfrc.amp{i}(:,:,1),36,19);


end