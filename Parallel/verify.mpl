list2maple := proc(ls :: list(nonnegint));
    local i;
    return add( seq( ls[i] * x^(i - 1), i = 1 .. nops(ls) ) );
end proc;


verifyFile := proc(in_name :: string, out_name :: string, skip :: posint := 0, $)
    local infile,outfile,i,a,b,c,p,res;
    
    p := 2^16 + 1;
    infile := fopen(in_name, READ);
    outfile := fopen(out_name, WRITE);

    try
        while not feof(infile) do
            a := list2maple(parse(readline(infile)));
            b := list2maple(parse(readline(infile)));
            c := list2maple(parse(readline(infile)));
            res := Expand(a*b) mod p;
            if res <> c then
                fprintf(outfile, "%a\n", [ seq( coeff(res, x, i), i=1..degree(res) ) ]);
            fi;

            # skip a bunch of lines
            for i to skip*3 do readline(infile); od;
        od;
    finally
        fclose(infile);
        fclose(outfile);
    end try;
end proc;


verifyFile("out", "actual");
