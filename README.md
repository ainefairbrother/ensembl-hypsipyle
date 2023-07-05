# GraphQL for Ensembl Variation

A [GraphQL trial](https://graphql.org/) for [Ensembl](https://www.ensembl.org) to reduce the need for RESTful web services.

This application is implemented with [Ariadne](https://ariadnegraphql.org/), a schema-first graphql framework for Python

GraphQL requires a schema (in /common) and implementation of resolver functions that know how to interpret specific parts of a GraphQL query. Resolvers are found in /resolver, and may also make use of "data loaders" to overcome inherent deficiencies in GraphQL implementations.


## Installation
Requires Python 3.8+.  

To install dependencies, run:

`pip install -r requirements.txt` for just the API.  Use this when deploying the service.

`pip install -r requirements-dev.txt` installs everything including dev dependencies like pytest, mypy etc.

## Running the API locally
Add MongoDB and GRPC configuration `./connections.conf`

The file follows the following template:
```
datafile=DATAFILE
```

This command will start the server:

```uvicorn --workers 1 --host=0.0.0.0 graphql_service.server:APP```


If you're developing in PyCharm, you will probably find it useful to create a run 
configuration so that you can use the debugger.  Create a run configuration that 
looks like this:


## Development

### Testing

TODO

### Linting

TODO

### Type checking

TODO

### Formatting

`black . --check --diff` previews the formatting.

`black .` applies the formatting in-place.

## Containerisation for dev

Build the image using `./Dockerfile.dev`:

`docker build -t $NAME:$VERSION -f ./Dockerfile.dev .`

Run a container with the image (`--publish` below is exposing the container's ports to the host network):

`docker container run --publish 0.0.0.0:80:80/tcp --publish 0.0.0.0:8000:8000/tcp -ti -v <ensembl-hysipile-dir>:/app $NAME:$VERSION`

The connection configuration is assumed to exist in the repo as the file `./connections.conf` and gets built into the Docker 
image. 


## Containerisation for prod
TODO