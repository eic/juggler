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
LONG_TAG=$(VERSION)-$(PUBLISH_TAG)

# help will output the help for each task
# thanks to https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help

help: ## This help.
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

.DEFAULT_GOAL := help

# ==========================================================================
#
build: ## build the image
	DOCKER_BUILDKIT=1 docker build -t $(APP_NAME):$(LONG_TAG) .

build-nc: ## Build the container without caching (from scratch)
	DOCKER_BUILDKIT=1 docker build --no-cache -t $(APP_NAME):$(LONG_TAG) .

# ==========================================================================
#
login: ## Auto login to AWS-ECR unsing aws-cli
	@docker login -u ${CI_REGISTRY_USER} -p ${CI_REGISTRY_PASSWORD} ${CI_REGISTRY}
	echo "Login COMPLETE"

# ==========================================================================
#

publish: login ## Publish a tagged container to ECR
	@echo 'publish $(PUBLISH_TAG) to $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME)'
	docker tag $(APP_NAME):$(LONG_TAG) $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(PUBLISH_TAG)
	docker push $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(PUBLISH_TAG)
	docker rmi  $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(PUBLISH_TAG) || true

publish-version: login ## Publish the `{version}` taged container to ECR
	@echo 'publish $(VERSION) to $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(VERSION)'
	docker tag $(APP_NAME):$(LONG_TAG) $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(VERSION)
	docker push $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(VERSION)
	docker rmi  $(REG_NAME)/$(GL_REG_GROUP)/$(APP_NAME):$(VERSION) || true

# ==========================================================================
# remove container from registry
unpublish: login
	@echo 'removing $(PUBLISH_TAG)'
	curl --request DELETE --header "PRIVATE-TOKEN: $(CI_JOB_TOKEN)" "$(REG_API_URL)/$(PUBLISH_TAG)"

# ==========================================================================
# cleanup docker registry on system used by runner
cleanup:
	@echo 'removing $(REG_NAME):$(LONG_TAG)'
	docker rmi $(REG_NAME):$(LONG_TAG) || true
