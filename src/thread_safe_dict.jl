using Base.Threads

mutable struct ThreadSafeDict{K,V}
    dict::Dict{K,V}
    lock::SpinLock

    ThreadSafeDict{K,V}() where {K,V} = new(Dict(), SpinLock())
end

# TODO: why am I not using wait() and notify()?
function Base.get!(default::Function, d::ThreadSafeDict{K,V}, k::K) where {K,V}
    #lock(d.lock)
    if !haskey(d.dict, k)
        d.dict[k] = default()
    end
    #unlock(d.lock)
    return d.dict[k]
end

function Base.push!(d::ThreadSafeDict{K,V}, kv::Pair{K,V}) where {K,V}
    #lock(d.lock)
    push!(d.dict, kv)
    #unlock(d.lock)
    return d
end

function Base.keys(d::ThreadSafeDict{K,V}) where {K,V}
    #lock(d.lock)
    k = keys(d.dict)
    #unlock(d.lock)
    return k
end

function Base.findmin(d::ThreadSafeDict{K,V}) where {K,V}
    #lock(d.lock)
    m = findmin(d.dict)
    #unlock(d.lock)
    return m
end

function Base.iterate(d::ThreadSafeDict{K,V}) where {K,V}
    #lock(d.lock)
    e = iterate(d.dict)
    #unlock(d.lock)
    return e
end

function Base.iterate(d::ThreadSafeDict{K,V}, i) where {K,V}
    #lock(d.lock)
    e = iterate(d.dict, i)
    #unlock(d.lock)
    return e
end

function Base.getindex(d::ThreadSafeDict{K,V}, k::K) where {K,V}
    #lock(d.lock)
    i = getindex(d.dict, k)
    #unlock(d.lock)
    return i
end