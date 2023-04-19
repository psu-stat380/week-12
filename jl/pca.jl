### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 79ef8d20-d9a7-11ed-06e0-35a08a52d6c9
using Pkg; Pkg.activate("/storage/home/suv87/work/courses/stat380/weeks/week-12/jl")

# ╔═╡ 30dc123d-a888-4778-83a9-3e819eef3356
begin
	using Plots, PlutoUI
	using LinearAlgebra, LazySets
	using Random, Distributions
	using Statistics, StatsBase
end

# ╔═╡ 4d1d9c70-5acb-49c3-9480-a4ef1e549eae
P(x, m, b=0.0) = (x[1] + m * x[2] - m * b, m * x[1] + m^2 * x[2] + b) ./ (1+m^2)

# ╔═╡ 1a04880f-4f33-471e-bb5a-2fead08ffad0
begin
	Σ = [1.0 -0.75; -0.75 1] .* 10
	μ = [0.0, 0.0]
	D = MvNormal(μ, Σ)
	X = rand(D, 300)
	Xn = Tuple.(eachcol(X))
end;

# ╔═╡ 08898403-e07f-44fd-8c20-ca8a84a37fe6
md"θ = $(@bind θ Slider(-pi + 1e-5 : 0.025 : pi - 1e-5, pi/4, true))"

# ╔═╡ 2215cc37-26c0-468e-8400-3e857b13b3ea
begin
	function pca(X)
		λ, v = eigen(cor(X, dims=2))
		round.(v[:, end], digits=2)
	end
	v1, v2= round.((cos(θ), sin(θ)), digits=2)
	m = tan(θ)
	Yn = map(x -> P(x, m), Xn)
	L = [LineSegment([x...], [y...]) for (x, y) in zip(Xn, Yn)]
end;

# ╔═╡ 0270e94b-f1b5-4a18-974b-b5658e8f56fa
begin
	gr(legend=:topright)
	
	######################################################
	plt1 = scatter(Xn, ma=0.2, label="", xlabel="x1", ylabel="x2", lim=(-10, 10), ratio=:equal, c=:dodgerblue)
		
	plt1 = plot(plt1, x -> tan(θ) * x, -10, 10, c=:black, lw=1, label="v = $((v1, v2))")
	# plt1 = plot(plt1, [(0.0, 0.0), (v1, v2) .* 4], arrow=true, lw=2, c=:firebrick1, label="v = $((v1, v2))")

	plt1 = plot(plt1, L, c=:red, la=0.5, ms=0)

	plt1 = scatter(plt1, Yn, label="", c=:red, ma=0.5, msw=0.5)	
	
	plt1 = scatter(plt1, Xn, ma=0.2, label="", c=:dodgerblue)
	
	
	######################################################
	V = map(x -> v1 * x[1] + v2 * x[2], Xn)
	
	plt2 = plot()
	
	plt2 = scatter(plt2, V, x -> 0, ylim=(-0.5, 0.5), xlim = (-20, 20), label="")
	title!("variance=$(round(var(V), digits=3))")
end;

# ╔═╡ 62021fb3-ce74-428e-9f4c-66fa498f3f32
# plot(plt1, size=(500, 500))
plot(plt1, plt2, size=(1000, 500))

# ╔═╡ ce9dccd1-4c09-40dc-b122-4858fa673ca4
pca(X)

# ╔═╡ 1e4e151b-3be4-42cf-9c95-bf413505eb00


# ╔═╡ cccbeaf1-fb23-4597-b778-52c4ea108527
html"<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>"

# ╔═╡ dbc4df55-3ec1-4c7a-9fa6-763d1f9923f4
md"# 3-d plot"

# ╔═╡ cef0bdc0-6d0f-46c0-9d9d-d8f0a39e1fe3
begin
	Σ3 = [1.0 0.75 0.5; 0.75 1.0 0.75; 0.5 0.75 1.0] .* 3
	μ3 = [0.0, 0.0, 0.0]
	D3 = MvNormal(μ3, Σ3)
	X3 = rand(D3, 500)
	Xn3 = Tuple.(eachcol(X3))
end;

# ╔═╡ fe06059e-7883-4abc-a47f-7d4128087dab
md"ϕ = $(@bind ϕ Slider(-pi + 1e-5 : 0.05 : pi - 1e-5, pi/4, true))",
md"ω = $(@bind ω Slider(-pi/2 + 1e-5 : 0.05 : pi/2 - 1e-5, pi/4, true))"

# ╔═╡ 0fa4d054-0897-48a6-9c99-8ca2b1cb08bf
begin
	u1, u2, u3 = round.( ( cos(ϕ) * cos(ω), sin(ϕ) * cos(ω), sin(ω) ), digits=2)
	H = Hyperplane([u1, u2, u3], 1.0)
end

# ╔═╡ 8a5d9774-8333-4466-b1a7-e8e704ddbf85
begin
	plotly()
	######################################################
	
	plt3 = scatter(Xn3, ms=2.0, label="", lim=(-10, 10), ratio=1, legend=:bottomright, ma=0.2, c=:black)
	
	# plt3 = plot(plt3, [(0.0, 0.0, 0.0), (u1, u2, u3) .* 5], c=:red, arrow=true, lw=5, ms=5, label="v = $((u1, u2, u3))")
	plt3 = plot(plt3, [(u1, u2, u3) .* -5 , (u1, u2, u3) .* 5], c=:red, arrow=true, lw=5, ms=5, label="v = $((u1, u2, u3))")

	######################################################
	V3 = map(x -> dot(x, (u1, u2, u3)), Xn3)
	
	plt4 = plot()

	plt4 = scatter(plt4, V3, x -> 0.0, ylim=(-1, 1), xlim=(-10, 10), label="variance = $(round(var(V3), digits=3))")
end;

# ╔═╡ f0610c39-8325-41ac-bb0a-d38fe8fb26e8
plot(plt3, size=(900, 500))
# plot(plt3, plt4, size=(900, 450))

# ╔═╡ f99b09eb-2abb-49db-9f8a-d94234a502cd
pca(X3)

# ╔═╡ Cell order:
# ╠═79ef8d20-d9a7-11ed-06e0-35a08a52d6c9
# ╟─30dc123d-a888-4778-83a9-3e819eef3356
# ╠═4d1d9c70-5acb-49c3-9480-a4ef1e549eae
# ╟─1a04880f-4f33-471e-bb5a-2fead08ffad0
# ╟─2215cc37-26c0-468e-8400-3e857b13b3ea
# ╠═0270e94b-f1b5-4a18-974b-b5658e8f56fa
# ╟─08898403-e07f-44fd-8c20-ca8a84a37fe6
# ╠═62021fb3-ce74-428e-9f4c-66fa498f3f32
# ╠═ce9dccd1-4c09-40dc-b122-4858fa673ca4
# ╠═1e4e151b-3be4-42cf-9c95-bf413505eb00
# ╟─cccbeaf1-fb23-4597-b778-52c4ea108527
# ╟─dbc4df55-3ec1-4c7a-9fa6-763d1f9923f4
# ╟─cef0bdc0-6d0f-46c0-9d9d-d8f0a39e1fe3
# ╟─0fa4d054-0897-48a6-9c99-8ca2b1cb08bf
# ╟─fe06059e-7883-4abc-a47f-7d4128087dab
# ╟─8a5d9774-8333-4466-b1a7-e8e704ddbf85
# ╠═f0610c39-8325-41ac-bb0a-d38fe8fb26e8
# ╠═f99b09eb-2abb-49db-9f8a-d94234a502cd
