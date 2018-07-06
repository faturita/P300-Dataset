% Export data

for s=subjectRange
    labelRange=SBJ(s).labelRange;
    csvwrite(sprintf('../BciSift/labelrange%d.csv',s),labelRange);
end
