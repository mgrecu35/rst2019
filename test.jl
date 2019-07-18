function test()
    ss=zeros(80)
    for i=1:100
        ss[1]+=1
    end
    println(ss[1])
end

test()
