"""
    Entity

A node of the IMO entity tree: `uuid`, `name`, `full_path` (e.g.
`"/satellite/LFT"`), `parent` UUID (`nothing` for root entities) and the UUIDs
of its `quantities`.
"""
struct Entity
	uuid::String
	name::String
	full_path::String
	parent::Union{Nothing, String}
	quantities::Vector{String}
end

"""
    Quantity

A quantity attached to an IMO entity: `uuid`, `name`, owning `entity` UUID,
`format_spec` UUID (or `nothing`) and the UUIDs of its `data_files`.
"""
struct Quantity
	uuid::String
	name::String
	entity::String
	format_spec::Union{Nothing, String}
	data_files::Vector{String}
end

"""
    DataFile

A version of a quantity in the IMO: `uuid`, `name`, `metadata` (the dictionary
the `DetectorInfo(imo, url)`-style constructors read), owning `quantity` UUID,
`upload_date`, `spec_version` and the `release_tags` it belongs to.
"""
struct DataFile
	uuid::String
	name::String
	metadata::Dict{String, Any}
	quantity::String
	upload_date::String
	spec_version::String
	release_tags::Vector{String}
end

"""
    Release

A tagged release of the IMO: `tag`, `release_date`, `comment` and the UUIDs of
the `data_files` it contains.
"""
struct Release
	tag::String
	release_date::String
	comment::String
	data_files::Set{String}
end

"""
    Imo(flatfile_location)

Instrument model database read from a local flat file, as `Imo` in litebird_sim
(through libinsdb). `flatfile_location` is either a directory containing
`schema.json` or the path of the schema file itself. Only uncompressed JSON is
supported: decompress `schema.json.gz` first.

Query objects with [`query`](@ref), or build parameters directly with
`DetectorInfo(imo, url)`, `FreqChannelInfo(imo, url)`, `InstrumentInfo(imo, url)`
and `SpinningScanningStrategy(imo, url)`.
"""
struct Imo
	location::String
	entities::Dict{String, Entity}
	path_to_entity::Dict{String, String}
	quantities::Dict{String, Quantity}
	data_files::Dict{String, DataFile}
	releases::Dict{String, Release}
end

function _schema_path(location::AbstractString)
	if isdir(location)
		path = joinpath(location, "schema.json")
		isfile(path) && return path
		if isfile(path * ".gz")
			throw(ArgumentError("$(path).gz is compressed; decompress it (e.g. `gunzip -k`) to read it with Falcons"))
		end
		throw(ArgumentError("no schema.json found in $location"))
	end
	isfile(location) || throw(ArgumentError("$location does not exist"))
	endswith(location, ".gz") && throw(ArgumentError("$location is compressed; decompress it to read it with Falcons"))
	return String(location)
end

_uuid_string(x) = lowercase(string(x))

function _walk_entities!(entities, objs, base_path, parent)
	for obj in objs
		uuid = _uuid_string(obj["uuid"])
		full_path = base_path * "/" * obj["name"]
		entities[uuid] = Entity(uuid, obj["name"], full_path, parent, String[])
		children = get(obj, "children", nothing)
		children === nothing || _walk_entities!(entities, children, full_path, uuid)
	end
	return entities
end

function Imo(flatfile_location::AbstractString)
	path = _schema_path(flatfile_location)
	schema = JSON.parsefile(path; dicttype = Dict{String, Any})

	entities = _walk_entities!(Dict{String, Entity}(), get(schema, "entities", []), "", nothing)
	path_to_entity = Dict(e.full_path => uuid for (uuid, e) in entities)

	quantities = Dict{String, Quantity}()
	for obj in get(schema, "quantities", [])
		uuid = _uuid_string(obj["uuid"])
		fmt = get(obj, "format_spec", nothing)
		quantities[uuid] = Quantity(uuid, get(obj, "name", ""), _uuid_string(obj["entity"]),
			fmt === nothing ? nothing : _uuid_string(fmt), String[])
		push!(entities[quantities[uuid].entity].quantities, uuid)
	end

	releases = Dict{String, Release}()
	tags_of = Dict{String, Vector{String}}()
	for obj in get(schema, "releases", [])
		files = Set(_uuid_string.(get(obj, "data_files", [])))
		releases[obj["tag"]] = Release(obj["tag"], string(get(obj, "release_date", "")),
			string(get(obj, "comment", get(obj, "comments", ""))), files)
		for uuid in files
			push!(get!(tags_of, uuid, String[]), obj["tag"])
		end
	end

	data_files = Dict{String, DataFile}()
	for obj in get(schema, "data_files", [])
		uuid = _uuid_string(obj["uuid"])
		metadata = something(get(obj, "metadata", nothing), Dict{String, Any}())
		data_files[uuid] = DataFile(uuid, get(obj, "name", ""), metadata, _uuid_string(obj["quantity"]),
			string(get(obj, "upload_date", "")), string(get(obj, "spec_version", "")), get(tags_of, uuid, String[]))
		push!(quantities[data_files[uuid].quantity].data_files, uuid)
	end
	return Imo(path, entities, path_to_entity, quantities, data_files, releases)
end

function Base.show(io::IO, imo::Imo)
	print(io, "Imo(\"", imo.location, "\", ", length(imo.data_files), " data files, releases ",
		join(sort(collect(keys(imo.releases))), ", "), ")")
end

const _UUID_RE = r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"

"The UUID at the end of a URL such as `/data_files/<uuid>`, or the string itself."
url_to_uuid(s::AbstractString) = _uuid_string(last(split(rstrip(s, '/'), '/')))

function _data_file(imo::Imo, uuid::AbstractString)
	key = _uuid_string(uuid)
	haskey(imo.data_files, key) || throw(ArgumentError("no data file with UUID $uuid in the IMO"))
	return imo.data_files[key]
end

"""
    query(imo, identifier) -> DataFile | Quantity | Entity

Look up an object of the IMO, as `Imo.query` in litebird_sim. `identifier` is

- a UUID string, or `"/data_files/<uuid>"`: a data file;
- `"/quantities/<uuid>"` or `"/entities/<uuid>"`: a quantity or an entity;
- a release path `"/releases/<tag>/<entity>/…/<quantity>"` (the `/releases`
  prefix is optional): the data file of that quantity included in the release.
  If the release contains several versions, the most recently uploaded one is
  returned.
"""
function query(imo::Imo, identifier::AbstractString)
	id = String(identifier)
	occursin(_UUID_RE, id) && return _data_file(imo, id)
	for (prefix, table) in (("/data_files/", imo.data_files), ("/quantities/", imo.quantities), ("/entities/", imo.entities))
		if startswith(id, prefix)
			uuid = url_to_uuid(id)
			haskey(table, uuid) || throw(ArgumentError("$identifier not found in the IMO"))
			return table[uuid]
		end
	end
	parts = filter(!isempty, split(startswith(id, "/releases/") ? id[length("/releases/")+1:end] : id, '/'))
	length(parts) >= 3 || throw(ArgumentError("malformed IMO path \"$identifier\"; expected /releases/<tag>/<entity>/…/<quantity>"))
	tag, entity_path, quantity_name = parts[1], "/" * join(parts[2:end-1], "/"), parts[end]
	haskey(imo.releases, tag) || throw(ArgumentError("release \"$tag\" not found; available: $(join(sort(collect(keys(imo.releases))), ", "))"))
	haskey(imo.path_to_entity, entity_path) || throw(ArgumentError("entity \"$entity_path\" not found in the IMO"))
	entity = imo.entities[imo.path_to_entity[entity_path]]
	qi = findfirst(q -> imo.quantities[q].name == quantity_name, entity.quantities)
	qi === nothing && throw(ArgumentError("entity \"$entity_path\" has no quantity \"$quantity_name\""))
	release = imo.releases[tag]
	candidates = [imo.data_files[u] for u in imo.quantities[entity.quantities[qi]].data_files if u in release.data_files]
	isempty(candidates) && throw(ArgumentError("quantity \"$entity_path/$quantity_name\" has no data file in release \"$tag\""))
	return candidates[argmax([c.upload_date for c in candidates])]
end

_metadata(imo::Imo, url::AbstractString) = query(imo, url).metadata
