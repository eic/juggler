# import config.
# You can change the default config with `make cnf="config_special.env" build`
cnf ?= config.env
include $(cnf)
# exports variables in config.env as environment variables
export $(shell sed 's/=.*//' $(cnf))

SHELL = bash
# grep the version from the mix file
VERSION=$(shell bash version.sh)
TAG_VERSION=$(VERSION)

# help will output the help for each task
# thanks to https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help

help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

.DEFAULT_GOAL := help

# ==========================================================================
#
build: ## build the image
	docker build -t $(APP_NAME) .

build-nc: ## Build the container without caching (from scratch)
	docker build --no-cache -t $(APP_NAME) .

# ==========================================================================
#
run: ## Run container on port configured in `config.env`
	docker run -i -t --rm --env-file=./config.env -p=$(PORT):$(PORT) --name="$(APP_NAME)" $(REPO)/$(APP_NAME):$(TAG_VERSION)

up: build run ## Run container on port configured in `config.env` (Alias to run)

stop: ## Stop and remove a running container
	docker stop $(APP_NAME); docker rm $(APP_NAME)

# ==========================================================================
#
tag: tag-latest tag-version #tag-version ## Generate container tags for the `{version}` and `latest` tags

tag-latest: ## Generate container `{version}` tag
	@echo 'create tag latest'
	#docker tag $(APP_NAME)   $(REPO)/$(APP_NAME):latest

tag-version: ## Generate container `latest` tag
	@echo 'creating tag $(APP_NAME):$(VERSION)'
	docker tag $(APP_NAME):latest $(APP_NAME):$(VERSION)


# ==========================================================================
#
login: ## Auto login to AWS-ECR unsing aws-cli
	@docker login -u ${CI_REGISTRY_USER} -p ${CI_REGISTRY_PASSWORD} ${CI_REGISTRY}
	echo "Login COMPLETE"

release-latest: build-nc login publish-latest ## Make a release by building and publishing the `{version}` ans `latest` tagged containers to ECR

release: build-nc publish clean-all## Make a release by building and publishing the `{version}` ans `latest` tagged containers to ECR

release-cached: build publish ## Make a release by building and publishing the `{version}` ans `latest` tagged containers to ECR

publish: login publish-latest publish-version #publish-version ## Publish the `{version}` ans `latest` tagged containers to ECR
	@echo "Publishing done"

publish-latest: ## Publish the `latest` taged container to ECR
	@echo 'publish latest to $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME)'
	docker tag $(APP_NAME):latest $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):latest
	docker push $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):latest
	docker rmi  $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):latest

publish-version: ## Publish the `{version}` taged container to ECR
	@echo 'publish latest to $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(VERSION)'
	docker tag $(APP_NAME):latest $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(VERSION)
	docker push $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(VERSION)
	docker rmi  $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(VERSION)

push: login ## push after login @echo 'push latest to $(REG_NAME)/$(REPO)/$(APP_NAME):latest'
	docker tag  $(APP_NAME):latest $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):latest
	docker push $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):latest
	docker rmi  $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):latest

# ==========================================================================
#

login-dockerhub: ## login to hub.docker.com
	@echo 'docker hub login :'
	docker login 

push-dockerhub: build #publish-version ## Publish the `{version}` ans `latest` tagged containers to ECR
	@echo '$(DH_ORG)/$(APP_NAME):$(TAG_VERSION)'
	docker push $(DH_ORG)/$(APP_NAME):latest
	#docker push $(DH_ORG)/$(APP_NAME):$(VERSION)

publish-dockerhub: build tag #publish-version ## Publish the `{version}` ans `latest` tagged containers to ECR
	@echo '$(DH_ORG)/$(APP_NAME):$(TAG_VERSION)'
	docker push $(DH_ORG)/$(APP_NAME):latest
	#docker push $(DH_ORG)/$(APP_NAME):$(VERSION)

# ==========================================================================
#
clean-tags:
	docker rmi  $(REPO)/$(APP_NAME):latest || true
	docker rmi  $(GL_GROUP)/$(APP_NAME):latest || true
	docker rmi  $(GL_REG_GROUP)/$(APP_NAME):latest || true
	docker rmi  $(REG_NAME)/$(REPO)/$(APP_NAME):latest   || true
	docker rmi  $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):latest || true

clean:
	docker rmi $(REPO_NAME)/$(APP_NAME) || true
	docker rmi $(REPO_NAME)/$(APP_NAME):$(TAG_VERSION) || true
	docker rmi $(REPO_NAME)/$(APP_NAME):$(VERSION) || true
	docker rmi $(REPO)/$(APP_NAME) || true
	docker rmi $(REPO)/$(APP_NAME):$(TAG_VERSION) || true
	docker rmi $(REPO)/$(APP_NAME):$(VERSION) || true
	docker rmi $(REG_NAME)/$(REPO)/$(APP_NAME):latest || true  
	docker rmi $(REG_NAME)/$(REPO)/$(APP_NAME):$(TAG_VERSION) || true  
	docker rmi $(REG_NAME)/$(REPO)/$(APP_NAME) || true  
	docker rmi $(DH_ORG)/$(APP_NAME) || true
	docker rmi $(DH_ORG)/$(APP_NAME):$(TAG_VERSION) || true
	docker rmi $(DH_ORG)/$(APP_NAME):$(VERSION) || true
	docker rmi $(GL_REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):latest || true  
	docker rmi $(GL_REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(VERSION) || true  
	docker rmi $(GL_REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME) || true  
	docker rmi $(GL_REG_GROUP)/$(APP_NAME) || true  
	docker rmi $(GL_REG_GROUP)/$(APP_NAME):$(VERSION) || true  
	#docker rmi $(APP_NAME) || true

clean-all: clean
	docker rmi $(REPO_NAME)/$(APP_NAME) || true

version: ## Output the current version
	@echo $(VERSION)

# ==========================================================================
#
info: ## Output the current version
	@echo 'VERSION      = $(VERSION)     '
	@echo 'REG_NAME     = $(REG_NAME)    '
	@echo 'APP_NAME     = $(APP_NAME)    '
	@echo 'REPO_NAME    = $(REPO_NAME)   '
	@echo 'DH_ORG       = $(DH_ORG)      '
	@echo 'GL_GROUP     = $(GL_GROUP)    '
	@echo 'GL_REG_GROUP = $(GL_REG_GROUP)'
	@echo 'GL_REG_NAME  = $(GL_REG_NAME) '
	@echo 'REPO_NAME    = $(REPO)        '
	@echo 'TAG_VERSION  = $(TAG_VERSION) '

ls: ## list all docker images 
	docker images 
