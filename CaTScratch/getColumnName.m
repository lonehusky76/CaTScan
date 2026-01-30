function colName = getColumnName(colNumber)
    colName = "";
    while colNumber > 0
        modulo = mod(colNumber - 1, 26);
        colName = char('A' + modulo) + colName;
        colNumber = floor((colNumber - modulo) / 26);
    end
end